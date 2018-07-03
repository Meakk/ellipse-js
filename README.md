# ellipse-js
A JavaScript library for ellipse fitting

## Example

This example computes equation of the best fitting ellipse of point set `(0,0), (1,0), (0,2), (2,1), (2,4)`

    var ellipse = new Ellipse();
    ellipse.setFromPoints([{x:0, y:0}, {x:1, y:0}, {x:0, y:2}, {x:2, y:1}, {x:2, y:4}]);
    ellipse.printEquation();

Prints `+0.743x^2 -0.557xy +0.371y^2 -0.743x -0.743y +0 = 0`

Each component can be accessed independently using `ellipse.a`, `ellipse.b`, `ellipse.c`, `ellipse.d`, `ellipse.e`.
