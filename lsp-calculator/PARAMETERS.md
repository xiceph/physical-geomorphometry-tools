# Land Surface Parameters

This page provides a list of Land Surface Parameters (LSPs) calculated by the `Land Surface Parameters Calculator` tool, along with detailed descriptions of each parameter's characteristics and equations.

All presented LSPs are local characteristics that can be defined (using differential geometry) in an infinitely small neighborhood of a given point. Most of them are functions of directional derivatives (changes) of elevation (zero-order characteristics) along the slope line and the contour line. LSPs are categorized into subsets based on the maximum order of derivatives used for their calculation. We denote the partial derivatives of the elevation z as follows $z_x = \frac{\partial z}{\partial x}$, $z_y = \frac{\partial z}{\partial y}$, $z_{xx} =\frac{\partial^2 z}{\partial x^2}$, $z_{xy} =\frac{\partial^2 z}{\partial x \partial y}$, $z_{yy} =\frac{\partial^2 z}{\partial x^2}$, $z_{xxx} =\frac{\partial^3 z}{\partial x^3}$, $z_{xxy} =\frac{\partial^3 z}{\partial x^2 \partial y}$, $z_{xyy} =\frac{\partial^3 z}{\partial x \partial y^2}$, $z_{yyy} =\frac{\partial^3 z}{\partial y^3}$.


## First-Order Parameters

Local first-order geomorphometric characteristics are those that can be defined using directional derivatives at a point, and their definition relies solely on the first partial derivatives of the elevation field $z_x$, $z_y$.

### Slope, Slope gradient – _S_

The slope gradient at a given point on a surface is defined as the angle between the tangent to the surface at that point (in the direction of the slope line) and the horizontal plane.

```math
S = \text{atan} \left( \sqrt{z_x^2 + z_y^2} \right)
```

#### Sine of Slope – sin _S_

The sine of the slope _S_ quantifies the potential gravitational energy per unit mass that can be converted into kinetic energy driving surface flows.

### Aspect – _A_

Aspect is the orientation of a slope (in the direction of descent), which indicates the direction of gravitational flows. It is the horizontal angle that depends on the orientation of the axes, traditionally measured clockwise from north. When using the easting-northing axes, the equation is (using the `atan2` function to determine the correct quadrant of the angle):

```math
A = \text{atan2}(z_y, z_x)
```

Due to the circular nature of aspect data, sine and cosine of the aspect are sometimes used for analysis (e.g., segmentation).

#### Sine of Aspect – sin _A_

The sine of the aspect _A_ represents the _y_-coordinate of the elevation gradient in the direction of decreasing elevation.

#### Cosine of Aspect – cos _A_

The cosine of the aspect _A_ represents the _x_-coordinate of the elevation gradient in the direction of decreasing elevation.

## Second-Order Parameters

Local second-order geomorphometric characteristics are defined by both the first ($z_x$, $z_y$) and second ($z_{xx}$, $z_{xy}$, $z_{yy}$) partial derivatives of the elevation field. Most of second-order characteristics describe the spatial changes in first-order characteristics. The vast majority of them describe various forms of surface curvature.

**Land surface curvatures** (LSC) represent a somewhat heterogeneous set of second-order geomorphometric characteristics, expressing curvature in different ways. A thorough analysis of this set of characteristics, in terms of their definition, calculation, and interpretation, was conducted in the work of Minár et al. (2020)[^1]. Here, we build on the latest version of this system.

![LSC](https://github.com/user-attachments/assets/25041a4d-9608-42dc-90a1-b679080c241e)



### Normal slope line (profile) curvature – (_k<sub>n</sub>_)_<sub>s</sub>_

Change of sine of slope (flow energy) in direction of slope line. Expression of energy of flow acceleration (change of flow energy to unit map distance). 

```math
(k_n)_s = -\frac{z_{xx} z_x^2 + 2z_{xy} z_x z_y + z_{yy} z_y^2}{(z_x^2 + z_y^2) \sqrt{(1 + z_x^2 + z_y^2)^3}}
```

### Normal contour (tangential) curvature – (_k<sub>n</sub>_)_<sub>c</sub>_

Divergence of slope lines multiplied by sine of slope. Expression of the energy of flow concentration.

```math
(k_n)_c = -\frac{z_{xx} z_y^2 + 2z_{xy} z_x z_y + z_{yy} z_x^2}{(z_x^2 + z_y^2) \sqrt{1 + z_x^2 + z_y^2}}
```

### Contour torsion – _𝜏<sub>c</sub>_

Change of slope in direction of contour line. Total measure of gravity discordant curvature contained in principal curvatures (_k<sub>max</sub>_ and _k<sub>min</sub>_) but not in gravity principal curvatures (_k<sub>n</sub>_)_<sub>s</sub>_ and (_k<sub>n</sub>_)_<sub>t</sub>_. 

```math
𝜏_c = \frac{z_x z_y (z_{xx} - z_{yy}) - z_{xy} (z_x^2 - z_y^2)}{(z_x^2 + z_y^2) \ (1 + z_x^2 + z_y^2)}
```

### Second slope line derivative - _z<sub>ss</sub>_

Second derivative of elevation, change of tangent of slope (map projection of the flow energy) in direction of slope line, expressing projected energy of flow acceleration.

```math
z_{ss} = \frac{(k_n)_s}{\cos S^3}
```

### Slope line torsion – _𝜏<sub>s</sub>_

Change of slope in direction of slope line. A projection of the energy of flow acceleration.

```math
𝜏_s = \frac{(k_n)_s}{\cos S}
```

### Second contour derivative - _z<sub>cc</sub>_

Second derivative of elevation. Projection of the energy of flow concentration to the map.

```math
z_{cc} = \frac{(k_n)_c}{\cos S}
```

### Projected contour (plan) curvature – (_k<sub>p</sub>_)_<sub>c</sub>_

Slope lines divergence → gravity flow divergence. 

```math
(k_p)_c = \frac{(k_n)_c}{\sin S}
```

### Projected slope line curvature (rotor) – (_k<sub>p</sub>_)_<sub>s</sub>_

Curvature of slope lines in map. Change of the aspect in direction of slope line. A measure of flow line twisting. 

```math
(k_p)_s = \frac{\text𝜏_c}{(\sin S \cdot \cos S)}
```

### Contour change of sine slope  – (sin _S_)_<sub>c</sub>_

General expression of gravity discordant (twisting) energy of flow (change of flow energy in contour direction).

```math
(\sin S)_c = \text𝜏_c \cdot \cos S
```

### Difference curvature - _k<sub>d</sub>_

Expression of summary energetic effect of profile and plan curvature on energy of gravity flow (local excess of flow energy).

```math
k_d = \frac{(k_n)_s - (k_n)_c}{2}
```

### Total accumulation curvature - _K<sub>a</sub>_

Product of gravity principal curvatures (_k<sub>n</sub>_)_<sub>s</sub>_ and (_k<sub>n</sub>_)_<sub>c</sub>_, expression of simplicity (linearity) at least in one of the gravity-principal directions.

```math
K_a = (k_n)_s \cdot (k_n)_c
```

### Total ring curvature - _K<sub>r</sub>_

Square of the contour torsion _𝜏<sub>c</sub>_, product of vertical and horizontal excess curvature _k<sub>ve</sub>_ and _k<sub>he</sub>_. Measure of flow line twisting.

```math
K_r = (k_n)_c \cdot (k_n)_s - k_{\text{min}} \cdot ((k_n)_c + (k_n)_s) + {k_\text{min}}^2
```

### Horizontal excess curvature - _k<sub>he</sub>_

Partial measure of gravity discordant curvature expressing difference between basic profile curvature (_k<sub>n</sub>_)_<sub>c</sub>_ and minimal curvature _k<sub>min</sub>_.

```math
k_{he} = (k_n)_c - k_{\text{min}}
```

### Vertical excess curvature - _k<sub>ve</sub>_

Partial measure of gravity discordant curvature expressing difference between basic plan curvature (_k<sub>n</sub>_)_<sub>s</sub>_ and minimal curvature _k<sub>min</sub>_.

```math
k_{ve} = (k_n)_s - k_{\text{min}}
```

### Maximal curvature - _k<sub>max</sub>_

Maximal curvature of normal section in given point determining first of the principal directions of surface deformation.

```math
k_{\text{max}} = - \frac{(1 + z_y^2) z_{xx} - 2 z_{xy} z_x z_y + (1 + z_x^2) z_{yy}}{2\sqrt{(1 + z_x^2 + z_y^2)^{3}}} + \sqrt{\left(\frac{(1 + z_y^2) z_{xx} - 2 z_{xy} z_x z_y + (1 + z_x^2) z_{yy}}{2\sqrt{(1 + z_x^2 + z_y^2)^{3}}}\right)^2 - \frac{z_{xx} z_{yy} - z_{xy}^2}{(1 + z_x^2 + z_y^2)^2}}
```

### Minimal curvature - _k<sub>min</sub>_

Minimal curvature of normal section in given point determining second of the principal directions of surface deformation.

```math
k_{\text{min}} = - \frac{(1 + z_y^2) z_{xx} - 2 z_{xy} z_x z_y + (1 + z_x^2) z_{yy}}{2 \sqrt{(1 + z_x^2 + z_y^2)^{3}}} - \sqrt{\left(\frac{(1 + z_y^2) z_{xx} - 2 z_{xy} z_x z_y + (1 + z_x^2) z_{yy}}{2\sqrt{(1 + z_x^2 + z_y^2)^{3}}}\right)^2 - \frac{z_{xx} z_{yy} - z_{xy}^2}{(1 + z_x^2 + z_y^2)^2}}
```

### Gaussian curvature - _K_

Product of principal curvatures _k<sub>max</sub>_ and _k<sub>min</sub>_, expression of simplicity (linearity) at least in one of principal directions. The only intrinsic curvature (detectable by a 2D observer) of the land surface.

```math
K = k_{\text{max}} \cdot k_{\text{min}}
```

### Elevation laplacian - _∇²z_

Double of mean of second derivatives of elevation. Map projection of diffusion energy of surface adequately describing the change of elevation with time.

```math
∇²z = \frac{(k_n)_c \cdot \cos S^2 + (k_n)_s}{\cos S^3}
```

### Unsphericity curvature  - _k<sub>u</sub>_

Difference of _k<sub>max</sub>_ and _k<sub>min</sub>_ expressing how far a surface is from a sphere (surface with minimal surface energy). Measure of general curvature anisotropy.

```math
k_u = \frac{k_{\text{max}} - k_{\text{min}}}{2}
```

### Mean curvature  - _k<sub>mean</sub>_

Mean of all normal curvatures. Measure of diffusion energy (potential local energy on the land surface applicable to diffusion).

```math
k_{\text{mean}} = \frac{(k_n)_s + (k_n)_c}{2} = \frac{k_{\text{max}} + k_{\text{min}}}{2}
```

### Casorati curvature  - _k<sub>c</sub>_

Measure of general curvedness of a surface, regardless of the sign of curvature (general nonlinearity) and integral expression of energies of all curvatures.

```math
k_c = \sqrt{\frac{k_{\text{max}}^2 + k_{\text{min}}^2}{2}}
```

## Third-Order Parameters

Local geomorphometric characteristics of the third order are defined by the first, second, and third partial derivatives of the elevation field. Since most second-order characteristics are referred to as curvatures, third-order characteristics, which describe their spatial variations, are termed changes of curvatures.

Due to the increased complexity of mathematical relationships in third-order parameters, we use the following substitutions for recurring elements:
```math
A = z_{xy}^2 - z_{xx} z_{yy}
```

```math
B = 2 z_x (z_{xxy} z_y + A) - z_{xxx} z_y^2 - z_{xyy} z_x^2 
```

```math
C = 2 z_y (z_{xyy} z_x + A) - z_{xxy} z_y^2 - z_{xxx} - z_{yyy} z_x^2 
```

```math
D = z_{xy} (z_{xx} + z_{yy})
```

```math
E = -2 (z_x (z_{xxy} z_y + z_{xy}^2 + z_{xx}^2)+ D z_y) - z_{xxx} z_x^2 - z_{xyy} z_y^2 
```

```math
F = -2 (z_y (z_{xyy} z_x + z_{xy}^2 + z_{yy}^2)+ D z_x) - z_{yyy} z_y^2 - z_{xxy} z_x^2 
```

```math
G = z_y^2 - z_x^2
```

```math
H = z_{xx} - z_{yy}
```

```math
I = z_{xxy} G+2 z_{xy} (z_{xy} z_y-z_{xx} z_x) + (z_{xx} z_y+z_{xy} z_x) H +z_x z_y (z_{xxx}-z_{xyy})
```

```math
J = z_{xyy} G+2 z_{xy} (z_{yy} z_y-z_{xy} z_x) + (z_{xy} z_y+z_{yy} z_x) H +z_x z_y (z_{xxy}-z_{yyy})
```

```math
K = 2 z_{xy} z_x z_y-z_{xx} z_y^2-z_{yy} z_x^2
```

```math
L = - z_x^2 z_{xx} - 2 z_{xy}z_xz_y - z_y^2-z_{yy}
```

```math
M = z_x^2+z_y^2
```

```math
N = z_{xx} z_x+z_{xy} z_y
```

```math
O = z_{yy} z_y+z_{xy} z_x
```

```math
P = 1+M
```

```math
R = z_x z_y H+z_{xy} G
```



### Contour change of normal contour curvature – (_k<sub>n</sub>_)_<sub>cc</sub>_

Expression of change of energy of flow concentration in contour direction.

```math
(k_n)_{cc} = \frac {z_y \left( BMP - K(2NP + MN)\right) - z_x \left( CMP - K(2OP + MO) \right)}{\sqrt{P^3} \sqrt{M^5}}
```

### Slope line change of normal contour curvature – (_k<sub>n</sub>_)_<sub>cs</sub>_

Expression of downstream change of energy of flow concentration.

```math
(k_n)_{cs} = \frac {z_x \left( K(2NP + MN) - BMP \right) - z_y \left( CMP - K(2OP + MO) \right)}{\sqrt{P^3} \sqrt{M^5}}
```

### Slope line change of normal slope line curvature – (_k<sub>n</sub>_)_<sub>ss</sub>_

Expression of downstream change of energy of flow acceleration (jerk or jolt).

```math
(k_n)_{ss} = \frac {P (z_x \left(2LN - EM \right) - z_y \left( FM - 2 LO \right)) + 3LM(z_y O + z_x N)}{\sqrt{P^5} \sqrt{M^5}}
```

[^1]: Minár, J., Evans, I. S., & Jenčo, M. (2020). A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews, 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
