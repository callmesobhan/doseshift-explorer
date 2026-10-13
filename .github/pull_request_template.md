## What and why

## Checklist
- [ ] `ruff check . && ruff format --check .` and `pytest` pass
- [ ] **No patient data** is included: no `*.dcm`, no image slices, nothing from `data/raw/` or `outputs/`
- [ ] New behaviour has a test (and a doc note if it changes a convention or a metric definition)
