# CalcShore Thermal Engine 2D

2D finite-difference thermal engine for mass concrete temperature control
plans (TCPs). Designed to produce output indistinguishable from
ConcreteWorks for contractor-facing reports.

## Status
Sprint 0 (2D port) in progress. See `docs/coding_passdown_v3.md`.

## Setup
    python -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    python -m pytest tests/ -v

## Milestones
- [ ] M0 — Grid builder + domain composition
- [ ] M1 — Pure conduction, analytical validation
- [ ] M2 — Hydration + variable props, centerline match vs. 1D v2
- [ ] M3 — Top BC, CW MIX-01 comparison
- [ ] M4 — Side BCs, ≤2°F RMS vs. CW
