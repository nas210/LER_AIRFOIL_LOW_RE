# PARSEC for Airfoils

## Introduction

In this repo, you will find an airfoil contour finder using modified PARSEC equations that account for pinned x and y values near the back of the airfoil. I wrote this code for a project investigating the effect leading-edge thickness has on LSB size and positioning at low reynolds while holding the forward thickness constant. 

## Background

PARSEC equations are a common way to describe the top and bottom curves of an airfoil. It allows for control of specific parameters under geometrical constraints. For more information about PARSEC equations, I recommended reading [Parametric Airfoils and Wings - Kimerius Aircraft.](https://www.kimerius.com/app/download/5782402238/Parametric+airfoils+and+wings.pdf)  

## Features

* Leading edge x and y coordinates: (`x_le`, `y_le`)
* Leading edge radius: (`r_le`)
* Trailing edge x and y coordinates: (`x_te`, `y_te`)
* Trailing edge angle: (`te_ang`)
* Crest x and y coordinates: (`x_c`, `y_c`)
* Slope at crest: (`dydx_c`)
* Curvature at crest: (`dy2dx2_c`)
* X pins: (`x_p1`,`x_p2`,`x_p3`)
* Y pins (y(x_p) = y_p): (`y_p1`,`y_p2`,`y_p3`)

## Dependencies
NumPy: Arrays
