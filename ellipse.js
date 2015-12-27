/*
The MIT License (MIT)

Copyright (c) 2015 Michael MIGLIORE

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

var Ellipse = (function() {

	// 3x3 matrix helpers
	function determinant(B) {
		return B[0][0] * B[1][1] * B[2][2]
		     + B[0][1] * B[1][2] * B[2][0]
		     + B[0][2] * B[1][0] * B[2][1]
		     - B[0][2] * B[1][1] * B[2][0]
		     - B[0][1] * B[1][0] * B[2][2]
		     - B[0][0] * B[1][2] * B[2][1];
	}

	function inverse(B) {
		var d = determinant(B);
		return [[(B[1][1] * B[2][2] - B[1][2] * B[2][1]) / d,
				 (B[0][2] * B[2][1] - B[0][1] * B[2][2]) / d,
				 (B[0][1] * B[1][2] - B[0][2] * B[1][1]) / d],
				[(B[1][2] * B[2][0] - B[1][0] * B[2][2]) / d,
				 (B[0][0] * B[2][2] - B[0][2] * B[2][0]) / d,
				 (B[0][2] * B[1][0] - B[0][0] * B[1][2]) / d],
				[(B[1][0] * B[2][1] - B[1][1] * B[2][0]) / d,
				 (B[0][1] * B[2][0] - B[0][0] * B[2][1]) / d,
				 (B[0][0] * B[1][1] - B[0][1] * B[1][0]) / d]];
	}

	function multiply(A, B) {
		return [[A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0],
				 A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1],
				 A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2]],
				[A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0],
				 A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1],
				 A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2]],
				[A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0],
				A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1],
				A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2]]];
	}

	function transpose(B) {
		return [[B[0][0], B[1][0], B[2][0]],
				[B[0][1], B[1][1], B[2][1]],
				[B[0][2], B[1][2], B[2][2]]];
	}

	function add(A, B) {
		return [[A[0][0] + B[0][0], A[0][1] + B[0][1], A[0][2] + B[0][2]],
				[A[1][0] + B[1][0], A[1][1] + B[1][1], A[1][2] + B[1][2]],
				[A[2][0] + B[2][0], A[2][1] + B[2][1], A[2][2] + B[2][2]]];
	}

	function trace(A) { return A[0][0] + A[1][1] + A[2][2]; }

	function scale(A, k) {
		return [[k * A[0][0], k * A[0][1], k * A[0][2]],
				[k * A[1][0], k * A[1][1], k * A[1][2]],
				[k * A[2][0], k * A[2][1], k * A[2][2]]];
	}
	
	function eigenvalues(A) {
		var q = trace(A) / 3;
		var K = add(A, [[-q, 0, 0],[0, -q, 0],[0, 0, -q]]);
		var p = Math.sqrt(trace(multiply(K,K))/6);
		var d = determinant(scale(K, 1 / p));

		var phi;
		if (d <= -2) {
			phi = Math.PI / 3;
		} else if (d >= 2) {
			phi = 0;
		} else {
			phi = Math.acos(d / 2) / 3;
		}
		
		return [q + 2 * p * Math.cos(phi),
				q + 2 * p * Math.cos(phi + (2 * Math.PI / 3)),
				q + 2 * p * Math.cos(phi + (4 * Math.PI / 3))];
	}
	
	function nullspace(G) {
		var k1 = -G[2][0]/G[2][2];
		var k2 = -G[2][1]/G[2][2];

		var y = -(G[1][0]+G[1][2]*k1)/(G[1][1]+G[1][2]*k2);
		var z = k1 + k2*y;
		var n = Math.sqrt(1+y*y+z*z);
		
		return [1/n, y/n, z/n];
	}

	var my = function(arg) {
		this.equation = {a:0, b:0, c:0, d:0, e:0, f:0, angle:0};
	
		this.setFromEquation = function(a, b, c, d, e, f) {
			this.equation.a = a;
			this.equation.b = b;
			this.equation.c = c;
			this.equation.d = d;
			this.equation.e = e;
			this.equation.f = f;
			this.equation.angle = 0;
		}
		
		this.setFromReducedEquation = function(a, c, d, e, f, angle) {
			this.equation.a = a;
			this.equation.b = 0;
			this.equation.c = c;
			this.equation.d = d;
			this.equation.e = e;
			this.equation.f = f;
			this.equation.angle = (angle === undefined)?0:angle;
		}
	
		this.setFromPoints = function(u){
			//compute sums
			var Sxxxx = u.reduce(function(p, c) { return p + c.x * c.x * c.x * c.x; }, 0);
			var Sxxxy = u.reduce(function(p, c) { return p + c.x * c.x * c.x * c.y; }, 0);
			var Sxxyy = u.reduce(function(p, c) { return p + c.x * c.x * c.y * c.y; }, 0);
			var Sxyyy = u.reduce(function(p, c) { return p + c.x * c.y * c.y * c.y; }, 0);
			var Syyyy = u.reduce(function(p, c) { return p + c.y * c.y * c.y * c.y; }, 0);
			var Sxxx  = u.reduce(function(p, c) { return p + c.x * c.x * c.x;       }, 0);
			var Sxxy  = u.reduce(function(p, c) { return p + c.x * c.x * c.y;       }, 0);
			var Sxyy  = u.reduce(function(p, c) { return p + c.x * c.y * c.y;       }, 0);
			var Syyy  = u.reduce(function(p, c) { return p + c.y * c.y * c.y;       }, 0);
			var Sxx   = u.reduce(function(p, c) { return p + c.x * c.x;             }, 0);
			var Sxy   = u.reduce(function(p, c) { return p + c.x * c.y;             }, 0);
			var Syy   = u.reduce(function(p, c) { return p + c.y * c.y;             }, 0);
			var Sx    = u.reduce(function(p, c) { return p + c.x;                   }, 0);
			var Sy    = u.reduce(function(p, c) { return p + c.y;                   }, 0);
			
			
			//construct matrices
			var S1 = [[Sxxxx, Sxxxy, Sxxyy],
					  [Sxxxy, Sxxyy, Sxyyy],
					  [Sxxyy, Sxyyy, Syyyy]];
			var S2 = [[Sxxx, Sxxy, Sxx],
					  [Sxxy, Sxyy, Sxy],
 					  [Sxyy, Syyy, Syy]];
			var S3 = [[Sxx, Sxy, Sx],
					  [Sxy, Syy, Sy],
  					  [Sx, Sy, u.length]];
			var S2T = transpose(S2);
			var iS3 = inverse(S3);
			var iC = [[0, 0, .5],
					  [0, -1, 0],
				 	  [.5, 0, 0]];

			var U = multiply(iS3, S2T);
			U = scale(U, -1);
			var A = multiply(iC, add(S1, multiply(S2, U)));
			
			var eig = eigenvalues(A);
			
			//get minimal positive eigenvalue
			var l = eig.reduce(function(p, c) {
				return (c > 0 && c < p) ? c : p;
			}, Infinity);
			
			//find eigenvector
			var ev = nullspace(add(A, [[-l, 0, 0],[0, -l, 0],[0, 0, -l]]));
			
			this.equation.a = ev[0];
			this.equation.b = ev[1];
			this.equation.c = ev[2];
			this.equation.d = U[0][0]*ev[0] + U[0][1]*ev[1] + U[0][2]*ev[2];
			this.equation.e = U[1][0]*ev[0] + U[1][1]*ev[1] + U[1][2]*ev[2];
			this.equation.f = U[2][0]*ev[0] + U[2][1]*ev[1] + U[2][2]*ev[2];
		}
		
		function printCoeff(x) {
			return (x<0?"-":"+") + Math.abs(Math.round(x*1000)/1000);
		}
		
		this.printEquation = function() {
			return printCoeff(this.equation.a) + "x^2 "
				 + printCoeff(this.equation.b) + "xy "
				 + printCoeff(this.equation.c) + "y^2 "
				 + printCoeff(this.equation.d) + "x "
				 + printCoeff(this.equation.e) + "y "
				 + printCoeff(this.equation.f) + " = 0";
		}
		
		this.convertToReducedEquation = function() {
			var eq = this.equation;
			var t = Math.atan(this.equation.b / (this.equation.c - this.equation.a))/2;
			var s = Math.sin(t);
			var c = Math.cos(t);
			var old_a = this.equation.a;
			var old_c = this.equation.c;
			var old_d = this.equation.d;
			var old_e = this.equation.e;
			this.equation.a = old_a*c*c - eq.b*c*s + old_c*s*s;
			this.equation.c = old_a*s*s + eq.b*c*s + old_c*c*c;
			this.equation.d = old_d*c - old_e*s;
			this.equation.e = old_d*s + old_e*c;
			this.equation.angle = t;
			this.equation.b = 0;
		}
		
		this.getAxisLength = function() {
			var eq = this.equation;
			if (Math.abs(eq.b) > 1e-9) this.convertToReducedEquation();
			var num = -4*eq.f*eq.a*eq.c + eq.c*eq.d*eq.d + eq.a*eq.e*eq.e;
			return [Math.sqrt(num/(4*eq.a*eq.c*eq.c)),
					Math.sqrt(num/(4*eq.a*eq.a*eq.c))];
		}
		
		this.getCenter = function() {
			var eq = this.equation;
			var denom = eq.b*eq.b - 4*eq.a*eq.c;
			return [(2*eq.c*eq.d - eq.b*eq.e)/denom,
					(2*eq.a*eq.e - eq.d*eq.b)/denom];
		}
	}
	return my;
})();
