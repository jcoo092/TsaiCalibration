module Datatypes
open MathNet.Numerics.LinearAlgebra

type CameraParams = {pixelSize: double; width: int; height: int; cx0: int; cy0: int}

type CalibrationPoint = {xw: int; yw: int; zw: int; cx: int; cy: int}

type Distortedcoords = {xd: double; yd: double}

type CalculatedParams = {R: Matrix<double>; T: Vector<double>; f: double; kappa1: double}

type CameraPoint = {xc: double; yc: double; zc: double}

type ImagePoint = {u: double; v: double; f: double}