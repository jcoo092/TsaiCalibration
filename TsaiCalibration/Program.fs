open System.IO
open MathNet.Numerics.LinearAlgebra

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

type CameraParams = {pixelSize: double; width: int; height: int}

let extractCameraParams cameraparamsfile = 
    match System.IO.Path.GetExtension(cameraparamsfile) with
    | ".csv" ->
        let rawLines = File.ReadAllLines(cameraparamsfile)
        {pixelSize = System.Double.Parse(rawLines.[0]); width = 640; height = 480}
    | _ -> failwith "invalid file type provided"

type CalibrationPoint = {xw: int; yw: int; zw: int; xd: int; yd: int}

let extractCalibrationPoints calibpointsfile = 
    match System.IO.Path.GetExtension(calibpointsfile) with
    | ".csv" ->
        let rawLines = File.ReadAllLines(calibpointsfile)
        let mutable splitLines = Array.map (fun (x : string) -> x.Split([|','|], System.StringSplitOptions.RemoveEmptyEntries)) rawLines
        if not (fst (System.Int32.TryParse(splitLines.[0].[0]))) then
            splitLines <- Array.skip 1 splitLines
        let nums = Array.map (fun x -> Array.map System.Int32.Parse x) splitLines
        Array.map (fun (x : int[]) -> {xw = x.[0]; yw = x.[1]; zw = x.[2]; xd = x.[3]; yd = x.[4]}) nums
    | _ -> failwith "invalid file type provided"

let openAndLoadFiles cameraparamsfilename calibrationpointsfilename = 
    let cameraparams = extractCameraParams cameraparamsfilename
    let calibrationpoints = extractCalibrationPoints calibrationpointsfilename
    (cameraparams, calibrationpoints)

[<EntryPoint>]
let main argv = 
    let cameraparams, calibrationpoints = openAndLoadFiles argv.[0] argv.[1]
    printfn "%A" argv
    0 // return an integer exit code
