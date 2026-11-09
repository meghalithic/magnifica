This repository is for the analysis of _Steginoporella magnifica_ images and traits.

# Images

The images are SEM images of _Steginoporella magnifica_ collected from New Zealand. These specimens are both modern and paleo (~2.3 Mya).

The images were collected by K. Voje, L.H. Liow, E. Di Martino, A. Piwoni-Piórewicz and others. The specimens were imaged by M. Ramsfjell.

Images are stored on a shared lab computer and will be made available with the publication of this project.

# Metadata

## Imaging metadata

The metadata file "[Imaged Steginoporella magnifica specimens.csv](https://github.com/megbalk/magnifica/blob/main/Imaged%20Steginoporella%20magnifica%20specimens.csv)" contains information about:
- date of image (Date) in MM/DD/YYY format
- a unique specimen number (specimenNR)
- number of pictures taken of the specimen (NR of pics)
- AV
- magnification
- formation
- comments

## Image protocol

There is also a file, "[Protocol Steginoporella magnifica.docx](https://github.com/megbalk/magnifica/blob/main/Protocol%20Steginoporella%20magnifica.docx)", containing information about the codes found in the metadata file.

## Location metadata

[Metadata for each formation](https://github.com/megbalk/stegino_metadata/blob/main/newMetadata/formations.csv) found in [Di Martino & Liow 2021](https://doi.org/10.1098/rspb.2020.2047).

# Traits

We extracted linear measurements from landmarks images of zooids.

## Landmarking

There are a total of 23 landmarks, numbered 0 to 22.

![landmarks](stegs_landmarks.png)

## Linear measurements

Linear measurements were extracted using the landmark coordinates.

![linear](stegs_linear_measurements.png)

*Traits*

- zooid length (zl; points 4 to 12)
- median process width (mpw.b; point 5 to 6)
- cryptocystal shelf width at midline (cw.m; points 10 to 11)
- cryptocystal shelf width at distal end (cw.d; point 8 to 7)
- cryptocystal side length (c.side; calculated as average between 1 to 7 and 9 to 8)
- orifice width at the midline (ow.m; points 19 to 0)
- orifice length (ol; calculated as hieght of triangle between points 4, 20, 21)
- orifice side length (o.side; calculated as average between 21 to 18 and 20 to 17) 

# Automation

A. Porto created [Steginator](https://github.com/agporto/Steginator), which links [DeepBryo](https://deepbryo.ngrok.io/), a tool developed by [DiMartino et al.](https://www.biorxiv.org/content/early/2022/11/17/2022.11.17.516938) that segments out zooids from a colony, and [ML-morph](https://github.com/agporto/ml-morph), a tool developed by [Porto & Voje](https://doi.org/10.1111/2041-210X.13373) to automatically place landmarks on images.

The output of the machine learning pipeline is a csv file of:
- id of the image (specimenNR_NR of pics_AV_magnification_backscatter)
- box_id, a unique identifier for the zooid in the image (box_top_box_left_box_width_box_height)
- box_top
- box_left
- box_width
- box_height
- Xn cooridnate, where n corresponds to the landmark number for landmarks 0-22
- Yn coordinate, where n corresponds to the landmark number for landmarks 0-22

# Data processing

In the code, "[fileNames.R](https://github.com/megbalk/magnifica/blob/main/fileNames.R)", reads in the image names and associated metadata file name and creates the dataset, "[imageList.csv](https://github.com/megbalk/magnifica/blob/main/imageList.csv)". This list retains information about the file structure, parses the image name into the same parts as the "Imaged Steginoporella magnifica specimens.csv", and notes the extension of the data item.

The code also checks for discrepencies being the imaging metadata, "Imaged Steginoporella magnifica specimens.csv", and the list of images.

# Exploratory analysis

The code, "[4-exploratoryAnalysis.R](https://github.com/megbalk/magnifica/blob/main/exploratoryAnalysis.R)", reads the file "[output.csv](https://github.com/megbalk/magnifica/blob/main/output.csv)", which is the output from the Steginator pipeline.

It extracts linear measurements using the landmark coordinates.

[metadata_check_BLEED.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/metadata_check_BLEED.R)  and [metadata_check_WABO_IV.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/metadata_check_WABO_IV.R) reconciles names and metadata of the images.

[1-filtering.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/1-filtering.R) ensures that the scale for the images are set and that all metadata is added.

# Scripts
[0-env.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/0-env.R) sets up the environment for downstream analyses.

[3-traits.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/3-traits.R) calculates linear measurements.

[5-gpMatrix.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/5-gpMatrix.R) calculates the **G** and **P** matrices.

[6-gpComparison.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/6-gpComparison.R) compares **G** and **P** and ∆z.

[temp.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/temp.R) calculates zooid length and temperature relationships.

[evo_analyses.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/evo_analyses.R) tests models of evolution.

[7-beta.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/7-beta.R) calculates the selection gradient and its affect on **G**.

[rarefaction_reimmanDistance.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/rarefaction_reimmanDistance.R) was an attempt to use Reimman Distances to test if **G** matrices were different or not based on sample size.

[plasticity.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/plasticity.R) was an attempt to calculate **E** and assess its influence on phenotypic change.

## Functions
[norm.vector.funct.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/norm.vector.funct.R) a function to normalize vectors.

[linearMeasurement.function.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/linearMeasurement.function.R) application of the Pythagorean Theorem to calculate linear distances.

[magnitudeFunction.R](https://github.com/meghalithic/magnifica/blob/main/Scripts/magnitudeFunction.R) a function to calculate the magnitude of beta.

