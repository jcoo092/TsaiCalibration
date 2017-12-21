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

type CalibrationPoint = {xw: int; yw: int; zw: int; xd: int; yd: int}

let extractCalibrationPoints calibpointsfile = 
    match System.IO.Path.GetExtension(calibpointsfile) with
    | ".csv" ->
        let rawLines = File.ReadAllLines(calibpointsfile)
        if rawLines.[0] System.String
        let splitLines = Array.map (fun (x : string) -> x.Split([|','|]) |> Array.map System.Int32.Parse) rawLines
        Array.map (fun (x : int[]) -> {xw = x.[0]; yw = x.[1]; zw = x.[2]; xd = x.[3]; yd = x.[4]}) splitLines
    | _ -> failwith "invalid file type provided"

let openAndLoadFiles cameraparamsfilename calibrationpointsfilename = 
    //let cameraparamsfile = File.Open(argv.[0], FileMode.Open, FileAccess.Read)
    //let calibrationpointsfile = File.Open(argv.[1], FileMode.Open, FileAccess.Read)
    let cameraparams = extractCameraParams cameraparamsfilename
    let calibrationpoints = extractCalibrationPoints calibrationpointsfilename
    (cameraparams, calibrationpoints)

[<EntryPoint>]
let main argv = 
    let cameraparams, calibrationpoints = openAndLoadFiles argv.[0] argv.[1]
    printfn "%A" argv
    0 // return an integer exit code
