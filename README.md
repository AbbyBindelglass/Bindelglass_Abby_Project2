# MATLAB Project 2

Estimate the arc length of the **upper boundary** of the Mandelbrot set with the following steps:
1) sampling boundary points with a **bisection** search in \(y\)
2) fitting a **15th order polynomial** to the **trimmed** curve segment
3) integrating \(\int \sqrt{1+(f'(x))^2}\,dx\) over the fitted interval

> Course: AMS 595 - "How Long is the Coast of Britain?"
> Author: Abby Bindelglass
> Due Date: 10/5
---
## Files
'Bindelglass_Abby_Project2.m' - main MATLAB script
'Bindelglass_Abby_Project2.pdf' - main MATLAB script converted to pdf
'Bindelglass_Abby_Project2_report.pdf' - report including work done, results, and discussion of results
---
## Requirements
- MATLAB
- No toolboxes needed
---
## How to Use
1. Open MATLAB
2. Run:
```matlab
Bindelglass_Abby_Project2.m
