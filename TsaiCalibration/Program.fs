open System.IO
open MathNet.Numerics.LinearAlgebra
open FSharp.Data

// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.

(*
Steps for Tsai Calibration:

1.  Read calibration points from file
2.  
*)

//let calculateBackProjectionErrorSum calibPoints = 
//    Array.map calculateBackProjectionError calibPoints |> Array.sum

//exception FileExtensionNotValid of string

type CameraParams = {pixelSize: double; width: int; height: int; cx0: int; cy0: int}

let extractCameraParams cameraparamsfile = 
    match System.IO.Path.GetExtension(cameraparamsfile) with
    | ".csv" ->
        let rawLines = File.ReadAllLines(cameraparamsfile)
        {pixelSize = System.Double.Parse(rawLines.[0]); width = 640; height = 480; cx0 = 320; cy0 = 240}
    | _ -> failwith "invalid file type provided"

type CalibrationPoint = {xw: int; yw: int; zw: int; cx: int; cy: int}

let extractCalibrationPoints calibpointsfile = 
    match System.IO.Path.GetExtension(calibpointsfile) with
    | ".csv" ->
        let rawLines = File.ReadAllLines(calibpointsfile)
        let mutable splitLines = Array.map (fun (x : string) -> x.Split([|','|], System.StringSplitOptions.RemoveEmptyEntries)) rawLines
        if not (fst (System.Int32.TryParse(splitLines.[0].[0]))) then
            splitLines <- Array.skip 1 splitLines
        let nums = Array.map (fun x -> Array.map System.Int32.Parse x) splitLines
        Array.map (fun (x : int[]) -> {xw = x.[0]; yw = x.[1]; zw = x.[2]; cx = x.[3]; cy = x.[4]}) nums
    | _ -> failwith "invalid file type provided"

let openAndLoadFiles cameraparamsfilename calibrationpointsfilename = 
    let cameraparams = extractCameraParams cameraparamsfilename
    let calibrationpoints = extractCalibrationPoints calibrationpointsfilename
    (cameraparams, calibrationpoints)

type Distortedcoords = {xd: double; yd: double}

let getdistortedcoords cameraparams calibrationpoints = 
    [|for calibpoint in calibrationpoints do
        yield {xd = cameraparams.pixelSize * double (calibpoint.cx - cameraparams.cx0); 
            yd = cameraparams.pixelSize * double (calibpoint.cy - cameraparams.cy0)}
    |]

let determineL cameraparams (calibrationpoints : CalibrationPoint[]) (distortedcoords : Distortedcoords[]) = 
    let n = calibrationpoints.Length
    let xds = Vector<double>.Build.DenseOfArray [| for i in 0 .. (n - 1) do yield distortedcoords.[n].xd |]
    //let yds = Vector<double>.Build.Dense n
    
    let mArr = [| for i in 0 .. (n - 1) do
                    let arr =  [| distortedcoords.[n].yd * double calibrationpoints.[n].xw;
                                    distortedcoords.[n].yd * double calibrationpoints.[n].yw;
                                    distortedcoords.[n].yd * double calibrationpoints.[n].zw;
                                    distortedcoords.[n].yd;
                                    -1.0 * distortedcoords.[n].xd * double calibrationpoints.[n].xw;
                                    -1.0 * distortedcoords.[n].xd * double calibrationpoints.[n].yw;
                                    -1.0 * distortedcoords.[n].xd * double calibrationpoints.[n].zw;
                    |]
                    yield arr |]

    let M = Matrix<double>.Build.DenseOfRowArrays(mArr)
    (M.TransposeThisAndMultiply M).PseudoInverse() * (M.Transpose() * xds)

let findDistanceFromCentre cameraParams calibrationPoint = 
    sqrt (double (calibrationPoint.cx - cameraParams.cx0) ** 2.0 + double (calibrationPoint.cy - cameraParams.cy0) ** 2.0)

let findSignOfTy cameraParams calibrationPoints (L : Vector<double>) absTy = 
    //let distancesFromCentre = Array.map (fun x -> sqrt (double (x.cx - cameraParams.cx0) ** 2.0 +
    //                                                double (x.cy - cameraParams.cy0) ** 2.0)) calibrationPoints |> Array.zip calibrationPoints
    let distancesFromCentre = Array.map (findDistanceFromCentre cameraParams) calibrationPoints |> Array.zip calibrationPoints
    let furthestPoint = Array.maxBy snd distancesFromCentre |> fst
    let r11 = L.[0] * absTy
    let r12 = L.[1] * absTy
    let r13 = L.[2] * absTy
    let r21 = L.[4] * absTy
    let r22 = L.[5] * absTy
    let r23 = L.[6] * absTy
    let tx = L.[3] * absTy
    let comp1 = r11 * double furthestPoint.xw + r12 * double furthestPoint.yw + r13 * double furthestPoint.zw + tx
    let comp2 = r21 * double furthestPoint.xw + r22 * double furthestPoint.yw + r23 * double furthestPoint.zw + absTy
    if System.Math.Sign comp1 = System.Math.Sign comp2 then
        absTy
    else
        -1.0 * absTy

let populateRandT (L : Vector<double>) sx (R : Matrix<double> ref) (T : Vector<double> ref) = 
    let ty = T.Value.[1]
    let ts = ty / sx
    R.Value.[0, 0] <- L.[0] * ts
    R.Value.[0, 1] <- L.[1] * ts
    R.Value.[0, 2] <- L.[2] * ts
    R.Value.[1, 0] <- L.[4] * ty
    R.Value.[1, 1] <- L.[5] * ty
    R.Value.[1, 2] <- L.[6] * ty
    T.Value.[0] <- L.[3] * ts

    // The below is supposed to represent taking the cross product of the first two rows of the rotation matrix

    let u = R.Value.[0, 0..2]
    let v = R.Value.[1, 0..2]
    R.Value.[2, 0] <- u.[1] * v.[2] - u.[2] * v.[1]
    R.Value.[2, 1] <- u.[2] * v.[0] - u.[0] * v.[2]
    R.Value.[2, 2] <- u.[0] * v.[1] - u.[1] * v.[0]

let calculateUy (R : Matrix<double> ref) (T : Vector<double> ref) cPoint = 
    R.Value.[1, 0] * double cPoint.xw + R.Value.[1, 1] * double cPoint.yw + R.Value.[1, 2] * double cPoint.zw + T.Value.[1]

let calculateUz (R : Matrix<double> ref) cPoint = 
    R.Value.[2, 0] * double cPoint.xw + R.Value.[2, 1] * double cPoint.yw + R.Value.[2, 2] * double cPoint.zw

let calculateFandTz calibrationpoints distortedcoords R T = 
    let Uy = Array.map (calculateUy R T) calibrationpoints
    let Uz = Array.map (calculateUz R) calibrationpoints
    let yds = [| for dc in distortedcoords do yield dc.yd |]
    let negatedYds = Array.map (fun x -> -1.0 * x) yds
    let M = Matrix<double>.Build.DenseOfRowArrays [Uy; negatedYds]
    let m = Vector<double>.Build.DenseOfArray (Array.map2 (fun x y -> x * y) Uz yds)
    let fandTz = (M.TransposeThisAndMultiply M).PseudoInverse() * (M.Transpose() * m)
    T.Value.[2] <- fandTz.[1]
    fandTz.[0]

[<EntryPoint>]
let main argv = 
    let cameraparams, calibrationpoints = openAndLoadFiles argv.[0] argv.[1]
    let distortedcoords = getdistortedcoords cameraparams calibrationpoints
    let L = determineL cameraparams calibrationpoints distortedcoords
    let absTy = (1.0 / (sqrt (L.[4] ** 2.0 + L.[5] ** 2.0 + L.[6] ** 2.0)))
    let sx = absTy * (sqrt (L.[0] ** 2.0 + L.[1] ** 2.0 + L.[2] ** 2.0))
    let R = ref (Matrix<double>.Build.Dense(3, 3, 0.0))
    let T = ref (Vector<double>.Build.Dense(3, 0.0))
    T.Value.[1] <- findSignOfTy cameraparams calibrationpoints L absTy // find final value of Ty
    populateRandT L sx R T
    let f = calculateFandTz calibrationpoints distortedcoords R T
    //printfn "%A" argv
    0 // return an integer exit code
