# STAGE1 — Diagnostic Temperature Trajectories

## Diagnostic Points

| Point | Width (m) | Depth (m) | Description |
| ----- | --------- | --------- | ----------- |
| side_near | 0.51 (idx 11) | 12.19 (idx 24) | Just inside concrete from side soil interface, mid-depth |
| bottom_near | 6.10 (idx 0) | 23.88 (idx 47) | Centerline, near bottom soil interface |
| deep_interior | 6.10 (idx 0) | 12.19 (idx 24) | Centerline, mid-depth |

## Figures

- `plots/trajectory_side_near.png`
- `plots/trajectory_bottom_near.png`
- `plots/trajectory_deep_interior.png`

## Symmetry Analysis

Max asymmetry `dT_pos + dT_neg` (should be 0 if linear):

- B↔D at Side-near (w≈0.51m, d≈12m): max asymmetry = 13.050°F
- B↔D at Bottom-near (centerline, d≈24m): max asymmetry = 13.050°F
- B↔D at Deep interior (centerline, d≈12m): max asymmetry = 13.068°F
- C↔E at Side-near (w≈0.51m, d≈12m): max asymmetry = 17.046°F
- C↔E at Bottom-near (centerline, d≈24m): max asymmetry = 17.046°F
- C↔E at Deep interior (centerline, d≈12m): max asymmetry = 17.064°F
- F↔H at Side-near (w≈0.51m, d≈12m): max asymmetry = 28.080°F
- F↔H at Bottom-near (centerline, d≈24m): max asymmetry = 28.080°F
- F↔H at Deep interior (centerline, d≈12m): max asymmetry = 28.098°F
- G↔I at Side-near (w≈0.51m, d≈12m): max asymmetry = 27.054°F
- G↔I at Bottom-near (centerline, d≈24m): max asymmetry = 27.054°F
- G↔I at Deep interior (centerline, d≈12m): max asymmetry = 27.090°F

## Saturation Temperature at t=168 hr

Expected: if system fully equilibrates, T_final → soil_temp_F.
Deviation indicates the thermal mass is still relaxing.

Run A (73°F/73°F) at Side-near (w≈0.51m, d≈12m): T_final=73.04°F
Run A (73°F/73°F) at Bottom-near (centerline, d≈24m): T_final=73.04°F
Run A (73°F/73°F) at Deep interior (centerline, d≈12m): T_final=73.13°F
Run B (73°F/60°F) at Side-near (w≈0.51m, d≈12m): T_final=63.97°F
Run B (73°F/60°F) at Bottom-near (centerline, d≈24m): T_final=63.97°F
Run B (73°F/60°F) at Deep interior (centerline, d≈12m): T_final=73.13°F
Run C (73°F/90°F) at Side-near (w≈0.51m, d≈12m): T_final=84.90°F
Run C (73°F/90°F) at Bottom-near (centerline, d≈24m): T_final=84.90°F
Run C (73°F/90°F) at Deep interior (centerline, d≈12m): T_final=73.13°F
Run D (60°F/73°F) at Side-near (w≈0.51m, d≈12m): T_final=69.10°F
Run D (60°F/73°F) at Bottom-near (centerline, d≈24m): T_final=69.10°F
Run D (60°F/73°F) at Deep interior (centerline, d≈12m): T_final=60.08°F
Run E (90°F/73°F) at Side-near (w≈0.51m, d≈12m): T_final=78.19°F
Run E (90°F/73°F) at Bottom-near (centerline, d≈24m): T_final=78.19°F
Run E (90°F/73°F) at Deep interior (centerline, d≈12m): T_final=90.16°F
Run F (73°F/45°F) at Side-near (w≈0.51m, d≈12m): T_final=53.51°F
Run F (73°F/45°F) at Bottom-near (centerline, d≈24m): T_final=53.51°F
Run F (73°F/45°F) at Deep interior (centerline, d≈12m): T_final=73.13°F
Run G (73°F/100°F) at Side-near (w≈0.51m, d≈12m): T_final=91.87°F
Run G (73°F/100°F) at Bottom-near (centerline, d≈24m): T_final=91.87°F
Run G (73°F/100°F) at Deep interior (centerline, d≈12m): T_final=73.13°F
Run H (45°F/73°F) at Side-near (w≈0.51m, d≈12m): T_final=64.54°F
Run H (45°F/73°F) at Bottom-near (centerline, d≈24m): T_final=64.54°F
Run H (45°F/73°F) at Deep interior (centerline, d≈12m): T_final=45.05°F
Run I (100°F/73°F) at Side-near (w≈0.51m, d≈12m): T_final=81.23°F
Run I (100°F/73°F) at Bottom-near (centerline, d≈24m): T_final=81.23°F
Run I (100°F/73°F) at Deep interior (centerline, d≈12m): T_final=100.18°F

## Linearity of ΔT Response

If cooling/heating rate scales linearly with |ΔT|, the ratio ΔT_actual/ΔT_input
should be constant across runs at the same diagnostic point.


### Side-near (w≈0.51m, d≈12m)
| Run | ΔT (soil-placement) | ΔT_actual at t=168 | Ratio |
| --- | ------------------- | ------------------ | ----- |
| B | -13°F | -9.072°F | 0.6978 |
| C | +17°F | +11.862°F | 0.6978 |
| D | +13°F | -3.942°F | -0.3032 |
| E | -17°F | +5.148°F | -0.3028 |
| F | -28°F | -19.530°F | 0.6975 |
| G | +27°F | +18.828°F | 0.6973 |
| H | +28°F | -8.496°F | -0.3034 |
| I | -27°F | +8.190°F | -0.3033 |

Linear fit slope: 0.1971 (1.0 = perfect linearity)
Max residual from linear fit: 14.0151°F

### Bottom-near (centerline, d≈24m)
| Run | ΔT (soil-placement) | ΔT_actual at t=168 | Ratio |
| --- | ------------------- | ------------------ | ----- |
| B | -13°F | -9.072°F | 0.6978 |
| C | +17°F | +11.862°F | 0.6978 |
| D | +13°F | -3.942°F | -0.3032 |
| E | -17°F | +5.148°F | -0.3028 |
| F | -28°F | -19.530°F | 0.6975 |
| G | +27°F | +18.828°F | 0.6973 |
| H | +28°F | -8.496°F | -0.3034 |
| I | -27°F | +8.190°F | -0.3033 |

Linear fit slope: 0.1971 (1.0 = perfect linearity)
Max residual from linear fit: 14.0151°F

### Deep interior (centerline, d≈12m)
| Run | ΔT (soil-placement) | ΔT_actual at t=168 | Ratio |
| --- | ------------------- | ------------------ | ----- |
| B | -13°F | +0.000°F | -0.0000 |
| C | +17°F | +0.000°F | 0.0000 |
| D | +13°F | -13.050°F | -1.0038 |
| E | -17°F | +17.028°F | -1.0016 |
| F | -28°F | +0.000°F | -0.0000 |
| G | +27°F | +0.000°F | 0.0000 |
| H | +28°F | -28.080°F | -1.0029 |
| I | -27°F | +27.054°F | -1.0020 |

Linear fit slope: -0.5012 (1.0 = perfect linearity)
Max residual from linear fit: 14.0457°F
