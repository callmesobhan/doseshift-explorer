"""Case discovery must never resolve ambiguity silently (scenarios from the independent review)."""

import copy
import shutil

import pydicom
import pytest
from pydicom.uid import generate_uid

from doseshift.cli import main
from doseshift.dicomrt import discover_cases, discover_with_report, write_case


@pytest.fixture
def folder(tmp_path, case):
    root = tmp_path / "case"
    write_case(case, root)
    return root


def _clone(src, dst, **changes):
    """Copy a DICOM file with a fresh SOP instance UID and some attribute changes."""
    ds = pydicom.dcmread(src)
    ds.SOPInstanceUID = generate_uid()
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    for k, v in changes.items():
        setattr(ds, k, v)
    ds.save_as(dst, enforce_file_format=True)
    return ds


def test_two_plans_on_one_structure_set_give_two_cases(folder):
    plan2 = _clone(folder / "RTPLAN.dcm", folder / "A_RTPLAN_boost.dcm", RTPlanLabel="BOOST")
    plan2.DoseReferenceSequence[0].TargetPrescriptionDose = 10.0
    plan2.save_as(folder / "A_RTPLAN_boost.dcm", enforce_file_format=True)
    dose = pydicom.dcmread(folder / "RTDOSE.dcm")
    d2 = _clone(folder / "RTDOSE.dcm", folder / "A_RTDOSE_boost.dcm")
    d2.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID = plan2.SOPInstanceUID
    d2.save_as(folder / "A_RTDOSE_boost.dcm", enforce_file_format=True)
    assert dose.SOPInstanceUID != d2.SOPInstanceUID
    cases = discover_cases(folder)
    assert len(cases) == 2
    assert sorted(c.prescription_gy for c in cases) == [10.0, 60.0]
    assert all(any("one case per plan" in n for n in c.notes) for c in cases)
    assert {c.rtdose.name for c in cases} == {"RTDOSE.dcm", "A_RTDOSE_boost.dcm"}


def test_plan_level_dose_beats_beam_doses_even_if_beams_sort_first(folder):
    for i in range(3):
        _clone(folder / "RTDOSE.dcm", folder / f"A_beam{i}.dcm", DoseSummationType="BEAM")
    (c,) = discover_cases(folder)
    assert c.rtdose.name == "RTDOSE.dcm"
    assert any("ignored 3 lower-priority RTDOSE" in n for n in c.notes)


def test_tied_doses_are_not_paired_and_say_why(folder):
    _clone(folder / "RTDOSE.dcm", folder / "RTDOSE_copy.dcm")  # same summation, same plan
    (c,) = discover_cases(folder)
    assert c.rtdose is None and not c.complete
    assert any("equally suitable RTDOSE" in n for n in c.notes)


def test_a_non_plan_level_only_dose_is_flagged(folder):
    ds = pydicom.dcmread(folder / "RTDOSE.dcm")
    ds.DoseSummationType = "BEAM"
    ds.save_as(folder / "RTDOSE.dcm", enforce_file_format=True)
    (c,) = discover_cases(folder)
    assert c.rtdose is not None and any(
        "is not plan-level (DoseSummationType BEAM)" in n for n in c.notes
    )


def test_structure_set_referencing_two_present_series_gives_two_cases(folder):
    ct = sorted((folder / "CT").glob("*.dcm"))
    other = folder / "OTHER"
    other.mkdir()
    new_uid = generate_uid()
    for f in ct:
        ds = pydicom.dcmread(f)
        ds.SeriesInstanceUID, ds.SeriesDescription = new_uid, "OTHER"
        ds.SOPInstanceUID = generate_uid()
        ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
        ds.save_as(other / f.name, enforce_file_format=True)
    rs = pydicom.dcmread(folder / "RTSTRUCT.dcm")
    extra = pydicom.dataset.Dataset()
    extra.SeriesInstanceUID = new_uid
    rs.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[
        0
    ].RTReferencedSeriesSequence.append(extra)
    rs.save_as(folder / "RTSTRUCT.dcm", enforce_file_format=True)
    cases = discover_cases(folder)
    assert len(cases) == 2 and {c.image_dir.name for c in cases} == {"CT", "OTHER"}
    assert all(any("one case per series" in n for n in c.notes) for c in cases)


def test_unplaceable_structure_sets_are_reported_not_dropped(folder):
    rs = pydicom.dcmread(folder / "RTSTRUCT.dcm")
    del rs.ReferencedFrameOfReferenceSequence
    rs.save_as(folder / "RTSTRUCT.dcm", enforce_file_format=True)
    cases, skipped = discover_with_report(folder)
    assert cases == [] and any("cannot be tied to an image" in s for s in skipped)
    with pytest.warns(UserWarning, match="skipped RTSTRUCT"):
        discover_cases(folder)


def test_cli_explains_why_nothing_was_found(folder, capsys):
    rs = pydicom.dcmread(folder / "RTSTRUCT.dcm")
    del rs.ReferencedFrameOfReferenceSequence
    rs.save_as(folder / "RTSTRUCT.dcm", enforce_file_format=True)
    assert main(["discover", str(folder)]) == 1
    out = capsys.readouterr().out
    assert "skipped RTSTRUCT" in out and "no analysable case" in out


def test_images_are_grouped_by_series_uid_not_by_which_file_sorts_first(folder):
    """A flat folder where another series' slice sorts first used to hide the structure set."""
    stray = _clone(sorted((folder / "CT").glob("*.dcm"))[0], folder / "CT" / "000_stray.dcm")
    stray.SeriesInstanceUID = generate_uid()
    stray.save_as(folder / "CT" / "000_stray.dcm", enforce_file_format=True)
    cases = discover_cases(folder)
    assert len(cases) == 1 and cases[0].image_dir.name == "CT"
    assert any("more than one image series" in n for n in cases[0].notes)


def test_multiple_target_doses_or_phases_withhold_the_prescription(folder):
    plan = pydicom.dcmread(folder / "RTPLAN.dcm")
    second = copy.deepcopy(plan.DoseReferenceSequence[0])
    second.DoseReferenceNumber, second.TargetPrescriptionDose = 2, 70.0
    plan.DoseReferenceSequence.append(second)
    plan.save_as(folder / "RTPLAN.dcm", enforce_file_format=True)
    (c,) = discover_cases(folder)
    assert c.prescription_gy is None and any("different target doses" in n for n in c.notes)
    plan = pydicom.dcmread(folder / "RTPLAN.dcm")
    del plan.DoseReferenceSequence[1]
    plan.FractionGroupSequence.append(copy.deepcopy(plan.FractionGroupSequence[0]))
    plan.save_as(folder / "RTPLAN.dcm", enforce_file_format=True)
    (c,) = discover_cases(folder)
    assert c.prescription_gy is None and any("multi-phase" in n for n in c.notes)


def test_structure_set_without_a_plan_says_a_plan_exists_elsewhere(folder):
    _clone(folder / "RTSTRUCT.dcm", folder / "RTSTRUCT_B.dcm")  # a recontour no plan references
    by_name = {c.rtstruct.name: c for c in discover_cases(folder)}
    b = by_name["RTSTRUCT_B.dcm"]
    assert b.rtplan is None and b.prescription_gy is None
    assert any("RTPLAN(s) present but none references this structure set" in n for n in b.notes)


def test_dose_in_another_frame_of_reference_is_flagged(folder):
    ds = pydicom.dcmread(folder / "RTDOSE.dcm")
    ds.FrameOfReferenceUID = generate_uid()
    ds.save_as(folder / "RTDOSE.dcm", enforce_file_format=True)
    (c,) = discover_cases(folder)
    assert any("FrameOfReference differs" in n for n in c.notes)


def test_symlinked_rt_files_are_found(folder, tmp_path):
    real = tmp_path / "elsewhere"
    real.mkdir()
    shutil.move(str(folder / "RTSTRUCT.dcm"), str(real / "RTSTRUCT.dcm"))
    (folder / "RTSTRUCT.dcm").symlink_to(real / "RTSTRUCT.dcm")
    assert len(discover_cases(folder)) == 1


# ---- residual cases from the third review ----
def test_a_dose_referencing_several_plans_is_matched_to_each(folder):
    dose = pydicom.dcmread(folder / "RTDOSE.dcm")
    other = pydicom.dataset.Dataset()
    other.ReferencedSOPClassUID = dose.ReferencedRTPlanSequence[0].ReferencedSOPClassUID
    other.ReferencedSOPInstanceUID = generate_uid()  # an unknown plan listed first
    dose.ReferencedRTPlanSequence.insert(0, other)
    dose.save_as(folder / "RTDOSE.dcm", enforce_file_format=True)
    (c,) = discover_cases(folder)
    assert c.rtdose is not None and not any("no RTDOSE references it" in n for n in c.notes)


def test_duplicate_files_of_one_object_give_one_case_and_a_note(folder):
    shutil.copy(folder / "RTSTRUCT.dcm", folder / "RTSTRUCT_copy.dcm")
    (c,) = discover_cases(folder)
    assert any(
        "duplicate file(s) of this structure set ignored: RTSTRUCT_copy.dcm" in n for n in c.notes
    )


def test_plan_level_but_effective_dose_gets_an_accurate_note(folder):
    ds = pydicom.dcmread(folder / "RTDOSE.dcm")
    ds.DoseType = "EFFECTIVE"
    ds.save_as(folder / "RTDOSE.dcm", enforce_file_format=True)
    (c,) = discover_cases(folder)
    assert any("plan-level but its DoseType is EFFECTIVE" in n for n in c.notes)
    assert not any("not plan-level" in n for n in c.notes)


def test_symlinked_series_directories_are_followed_and_cycles_do_not_hang(folder, tmp_path):
    real = tmp_path / "elsewhere" / "CT"
    real.parent.mkdir()
    shutil.move(str(folder / "CT"), str(real))
    (folder / "CT").symlink_to(real, target_is_directory=True)
    (folder / "loop").symlink_to(folder, target_is_directory=True)  # a cycle
    cases = discover_cases(folder)
    assert len(cases) == 1 and cases[0].complete
