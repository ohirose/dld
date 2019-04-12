## Simple instruction for running this software

This is the directory for running a demo program for a face dataset.
If you are a MATLAB user in Windows 10, you can run the demo program by simply typing

- trainingFace (+Enter)
- runDemoFace  (+Enter)

after changing your working directory to this directory in the MATLAB command window.
The training will be executed in a leave-one-person-out manner, i.e., leave six faces out.

This demo program registers the average shape and a face included in the face dataset.
If you would like to try another face shape, modify the variables `HUMAN_ID`, `FACE_ID`
and `SEX`  in `runDemoFace.m` into appropriate ones, each of which can be checked at
the 'Data' directory.

For more information such as optional parameters, see `README.md` above this directory,
please.
