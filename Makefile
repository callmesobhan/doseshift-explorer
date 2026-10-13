.PHONY: install test lint format demo figures app clean data real

install:
	pip install -e ".[all]"

test:
	pytest

lint:
	ruff check . && ruff format --check .

format:
	ruff format . && ruff check . --fix

demo:
	doseshift demo --out outputs/demo

# refresh the README gallery from a fresh demo run
figures: demo
	cp outputs/demo/grid_audit/geometry.png docs/img/geometry.png
	cp outputs/demo/sweep/sweep.png docs/img/sweep.png
	cp outputs/demo/grid_audit/audit_equivalence.png docs/img/audit_equivalence.png
	cp outputs/demo/grid_audit/dvh_scenarios.png docs/img/dvh_scenarios.png
	cp outputs/demo/monte_carlo/mc_dvh_band.png docs/img/mc_dvh_band.png
	cp outputs/demo/margin_scan/margin_scan.png docs/img/margin_scan.png

# ~0.6 GB; run the scripts without --yes first to see sizes and free disk
data:
	python scripts/download_tcia.py Vestibular-Schwannoma-SEG data/raw/vs-seg VS-SEG-135 VS-SEG-154 VS-SEG-188 --yes
	python scripts/download_tcia.py Pancreatic-CT-CBCT-SEG data/raw/pancreas-ct-cbct-seg Pancreas-CT-CB_014 Pancreas-CT-CB_017 --yes

real:
	for c in experiments/real/*.yaml; do doseshift run $$c --out outputs/real/$$(basename $$c .yaml); done

app:
	streamlit run app/streamlit_app.py

clean:
	rm -rf outputs .pytest_cache .ruff_cache
