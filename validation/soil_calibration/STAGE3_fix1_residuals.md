# STAGE3 Fix 1 — Residual Summary

**Fix applied:** soil_temp_F plumbed; is_submerged=False (concrete sides still air)  
**Comparison time:** t=168 hr  
**Method:** bilinear interpolation of engine (21×13) onto CW (49×13) grid.

## Residual Table

| Run | Placement/Soil | max\|R\| (°F) | mean\|R\| (°F) | RMS R (°F) | Location of max |
| --- | -------------- | -------------- | --------------- | ------------ | --------------- |
| A | 73/73 | 1.655 | 0.406 | 0.617 | depth=2.03m, w=0.00m |
| B | 73/60 | 11.436 | 2.764 | 4.511 | depth=0.00m, w=0.00m |
| C | 73/90 | 18.647 | 4.262 | 7.064 | depth=2.03m, w=0.00m |
| D | 60/73 | 14.581 | 3.306 | 5.502 | depth=2.54m, w=0.00m |
| E | 90/73 | 15.292 | 3.677 | 6.021 | depth=0.00m, w=0.00m |
| F | 73/45 | 26.448 | 6.250 | 10.292 | depth=0.00m, w=0.00m |
| G | 73/100 | 28.655 | 6.578 | 10.912 | depth=2.03m, w=0.00m |
| H | 45/73 | 29.498 | 6.749 | 11.246 | depth=4.06m, w=0.00m |
| I | 100/73 | 25.209 | 5.988 | 9.853 | depth=13.21m, w=0.00m |
