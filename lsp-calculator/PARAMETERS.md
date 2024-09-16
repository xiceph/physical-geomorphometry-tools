# Land Surface Parameters

This page provides a list of Land Surface Parameters (LSPs) calculated by the `Land Surface Parameters Calculator` tool, along with detailed descriptions of each parameter's characteristics and equations.

All presented LSPs are local characteristics that can be defined (using differential geometry) in an infinitely small neighborhood of a given point. They represent directional derivatives (changes) of elevation (zero-order characteristics) along the slope line and the contour line. LSPs are categorized into subsets based on the maximum order of derivatives used for their calculation.


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

Local second-order geomorphometric characteristics are defined by both the first ($z_x$, $z_y$) and second ($z_{xx}$, $z_{xy}$, $z_{yy}$) partial derivatives of the elevation field. Second-order characteristics describe the spatial changes in first-order characteristics. The vast majority of them describe various forms of surface curvature.

**Land surface curvatures** (LSC) represent a somewhat heterogeneous set of second-order geomorphometric characteristics, expressing curvature in different ways. A thorough analysis of this set of characteristics, in terms of their definition, calculation, and interpretation, was conducted in the work of Minár et al. (2020)[^1]. Here, we build on the latest version of this system.

![LSC](https://github.com/user-attachments/assets/3e47b184-37a7-46a6-86a5-7cc10621028e)


### Normal slope line (profile) curvature – (_k<sub>n</sub>_)_<sub>s</sub>_

### Normal contour (tangential) curvature – (_k<sub>n</sub>_)_<sub>c</sub>_

### Contour (geodesic) torsion – _𝜏<sub>c</sub>_


### Second slope line derivative

### Slope line torsion

### Second contour derivative

### Projected contour curvature

### Projected slope line curvature

### Contour change of sine slope


### Difference curvature

### Total accumulation curvature

### Total ring curvature

### Horizontal excess curvature

### Vertical excess curvature


### Maximal curvature

### Minimal curvature


### Gaussian curvature

### Elevation laplacian

### Unsphericity curvature

### Mean curvature

### Casorati curvature

## Third-Order Parameters

### Contour change of normal contour curvature

### Slope line change of normal contour curvature

### Slope line change of normal slope line curvature

[^1]: Minár, J., Evans, I. S., & Jenčo, M. (2020). A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews, 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
