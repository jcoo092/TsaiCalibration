module Projections
open MathNet.Numerics.LinearAlgebra
open Datatypes
open MathNet.Numerics.Optimization

// Forward projection - world to image plane

let worldToCameraTransform calculatedparams calibrationpoint = 
    let worldVec = Vector<double>.Build.DenseOfArray([|double calibrationpoint.xw; double calibrationpoint.yw; double calibrationpoint.zw|])
    calculatedparams.R.Multiply(worldVec) + calculatedparams.T

let cameraToImagePlaneTransform (calculatedparams : CalculatedParams) (camerapoint : Vector<double>) = 
    let u = calculatedparams.f * (camerapoint.[0] / camerapoint.[2])
    let v = calculatedparams.f * (camerapoint.[1] / camerapoint.[2])
    {u = u; v = v; f = calculatedparams.f}

let projectWorldToImage calculatedparams calibrationpoint = 
    worldToCameraTransform calculatedparams calibrationpoint |> cameraToImagePlaneTransform calculatedparams

let projectArrayWorldToImage calculatedparams calibrationpoints = 
    Array.map (projectWorldToImage calculatedparams) calibrationpoints

let calculateSADProjectionError (calibrationpoint, imagepoint) = 
    System.Math.Abs(double calibrationpoint.cx - imagepoint.u) + System.Math.Abs(double calibrationpoint.cy - imagepoint.v)

let calculateSSDProjectionError (calibrationpoint, imagepoint) = 
    System.Math.Pow(double calibrationpoint.cx - imagepoint.u, 2.0) + System.Math.Pow(double calibrationpoint.cy - imagepoint.v, 2.0)

let calculateProjectionErrors errorFunction calculatedparams calibrationpoints = 
    projectArrayWorldToImage calculatedparams calibrationpoints |> Array.zip calibrationpoints |> Array.map errorFunction

let calculateProjectionErrorsSAD calculatedparams calibrationpoints = 
    calculateProjectionErrors calculateSADProjectionError calculatedparams calibrationpoints

let calculateProjectionErrorsSSD calculatedparams calibrationpoints = 
    calculateProjectionErrors calculateSSDProjectionError calculatedparams calibrationpoints

// back projection - image plane to world



// Optimisation - just sticking this here for the time being
// http://www.imagingshop.com/linear-and-nonlinear-least-squares-with-math-net/ seems to have a reasonably good write up of
// optimisation with Math.Net

// Numerical differentiation section:

// numerical estimation of differentation, as straight out of Hughes' 1987 paper 'Why functional programming matters'
// (see https://doi.org/10.1093/comjnl/32.2.98).  See also https://www.youtube.com/watch?v=vGVJYoKIzjU where Hughes
// re-summarises the paper in a modern retrospective.

// Function to create a stream of successive evaluations of an input function
// Using the Seq type here means that the stream is lazily evaluated
let rec repeat f a = seq {yield a; yield! repeat f (f a)}

// Function to find the value in a sequence where the difference between that value and the previous value is smaller than a given tolerance
let rec within eps values = 
    let a = Seq.head values
    let b = Seq.head (Seq.tail values)
    if abs(a - b) <= eps then
        b
    else
        Seq.tail values |> within eps

// Function to estimate the first derivative of a given single-variable numerical function, at any given point for which the function is defined.
// Note that h should very small here (but be careful of precision issues when working with floating point numbers)
let easydiff f x h = ((f (x+h)) - f x) / h

// Helper function used in the succession of approximations
let halve x = x / 2.0

// Helper function that returns a sequence of better approximations of the differential
let differentiate h0 f x = repeat halve h0 |> Seq.map (easydiff f x) 

// A function used to help eliminate the error term that is a part of the differentiation estimation used above
// This function is used in between the differentiate function which generates the normal estimate sequence, and the
// within function, which checks if there is a sufficiently small change between estimates
// Note that this function is specific to the differentiation technique used
let rec elimerror n values = 
    let a = Seq.item 0 values
    let b = Seq.item 1 values
    seq {yield ((b * System.Math.Pow(2.0, n) - a) / (System.Math.Pow(2.0, n) - 1.0)); yield! elimerror n (Seq.tail values)}

let rec elimall n s = seq {yield s; yield! elimall (n + 1.0) (elimerror n s)}

// Takes the first term from the above sequence of sequences, as required
// In other words, it gets more and more precise values as we go on
let super s = Seq.map Seq.head (elimall 1.0 s)

// A derivation function, which should be able to calculate the differential at a given point
// (so long as a derivative is mathematically defined for the given input to the derived function)
// eps is the tolerance at which the function stops, when the difference between steps is too small
// h0 is the 'initial guess', which can probably be generally set to 1.0 to be honest
// f is the function that we seek the differential of
// x is the input to f that we want the numerical differential for
let deriveF eps h0 f x = differentiate h0 f x |> super |> within eps