%% About defineImage.mlx
% This file defines the MATLAB interface to the library |Image|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for &lt;SHAPE&gt;, &lt;DIRECTION&gt;, etc. For more
% information, see <matlab:helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') Define MATLAB Interface for C++ Library>.



%% Setup
% Do not edit this setup section.
function libDef = defineImage()
libDef = clibgen.LibraryDefinition("ImageData.xml");
%% OutputFolder and Libraries 
libDef.OutputFolder = "E:\code\functionalimaginganalysis\Image";
libDef.Libraries = "E:\code\functionalimaginganalysis\Image\Image.dll";

%% C++ class |Vec3<int>| with MATLAB name |clib.Image.Vec3_int_| 
Vec3_int_Definition = addClass(libDef, "Vec3<int>", "MATLABName", "clib.Image.Vec3_int_", ...
    "Description", "clib.Image.Vec3_int_    Representation of C++ class Vec3<int>."); % Modify help description values as needed.

%% C++ class constructor for C++ class |Vec3<int>| 
% C++ Signature: Vec3<int>::Vec3(Vec3<int> const & input1)
Vec3_int_Constructor1Definition = addConstructor(Vec3_int_Definition, ...
    "Vec3<int>::Vec3(Vec3<int> const & input1)", ...
    "Description", "clib.Image.Vec3_int_ Constructor of C++ class Vec3<int>."); % Modify help description values as needed.
defineArgument(Vec3_int_Constructor1Definition, "input1", "clib.Image.Vec3_int_", "input");
validate(Vec3_int_Constructor1Definition);

%% C++ class constructor for C++ class |Vec3<int>| 
% C++ Signature: Vec3<int>::Vec3()
Vec3_int_Constructor2Definition = addConstructor(Vec3_int_Definition, ...
    "Vec3<int>::Vec3()", ...
    "Description", "clib.Image.Vec3_int_ Constructor of C++ class Vec3<int>."); % Modify help description values as needed.
validate(Vec3_int_Constructor2Definition);

%% C++ class public data member |x| for C++ class |Vec3<int>| 
% C++ Signature: int Vec3<int>::x
addProperty(Vec3_int_Definition, "x", "int32", ...
    "Description", "int32    Data member of C++ class Vec3<int>."); % Modify help description values as needed.

%% C++ class public data member |y| for C++ class |Vec3<int>| 
% C++ Signature: int Vec3<int>::y
addProperty(Vec3_int_Definition, "y", "int32", ...
    "Description", "int32    Data member of C++ class Vec3<int>."); % Modify help description values as needed.

%% C++ class public data member |z| for C++ class |Vec3<int>| 
% C++ Signature: int Vec3<int>::z
addProperty(Vec3_int_Definition, "z", "int32", ...
    "Description", "int32    Data member of C++ class Vec3<int>."); % Modify help description values as needed.

%% C++ class |Vec3<double>| with MATLAB name |clib.Image.Vec3_double_| 
Vec3_double_Definition = addClass(libDef, "Vec3<double>", "MATLABName", "clib.Image.Vec3_double_", ...
    "Description", "clib.Image.Vec3_double_    Representation of C++ class Vec3<double>."); % Modify help description values as needed.

%% C++ class constructor for C++ class |Vec3<double>| 
% C++ Signature: Vec3<double>::Vec3(Vec3<double> const & input1)
Vec3_double_Constructor1Definition = addConstructor(Vec3_double_Definition, ...
    "Vec3<double>::Vec3(Vec3<double> const & input1)", ...
    "Description", "clib.Image.Vec3_double_ Constructor of C++ class Vec3<double>."); % Modify help description values as needed.
defineArgument(Vec3_double_Constructor1Definition, "input1", "clib.Image.Vec3_double_", "input");
validate(Vec3_double_Constructor1Definition);

%% C++ class constructor for C++ class |Vec3<double>| 
% C++ Signature: Vec3<double>::Vec3()
Vec3_double_Constructor2Definition = addConstructor(Vec3_double_Definition, ...
    "Vec3<double>::Vec3()", ...
    "Description", "clib.Image.Vec3_double_ Constructor of C++ class Vec3<double>."); % Modify help description values as needed.
validate(Vec3_double_Constructor2Definition);

%% C++ class public data member |x| for C++ class |Vec3<double>| 
% C++ Signature: double Vec3<double>::x
addProperty(Vec3_double_Definition, "x", "double", ...
    "Description", "double    Data member of C++ class Vec3<double>."); % Modify help description values as needed.

%% C++ class public data member |y| for C++ class |Vec3<double>| 
% C++ Signature: double Vec3<double>::y
addProperty(Vec3_double_Definition, "y", "double", ...
    "Description", "double    Data member of C++ class Vec3<double>."); % Modify help description values as needed.

%% C++ class public data member |z| for C++ class |Vec3<double>| 
% C++ Signature: double Vec3<double>::z
addProperty(Vec3_double_Definition, "z", "double", ...
    "Description", "double    Data member of C++ class Vec3<double>."); % Modify help description values as needed.

%% C++ class |Image<short>| with MATLAB name |clib.Image.Image_short_| 
Image_short_Definition = addClass(libDef, "Image<short>", "MATLABName", "clib.Image.Image_short_", ...
    "Description", "clib.Image.Image_short_    Representation of C++ class Image<short>."); % Modify help description values as needed.

%% C++ class method |shape| for C++ class |Image<short>| 
% C++ Signature: Vec3<int> Image<short>::shape()
shapeDefinition = addMethod(Image_short_Definition, ...
    "Vec3<int> Image<short>::shape()", ...
    "MATLABName", "shape", ...
    "Description", "shape Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(shapeDefinition, "RetVal", "clib.Image.Vec3_int_");
validate(shapeDefinition);

%% C++ class method |nx| for C++ class |Image<short>| 
% C++ Signature: int Image<short>::nx()
nxDefinition = addMethod(Image_short_Definition, ...
    "int Image<short>::nx()", ...
    "MATLABName", "nx", ...
    "Description", "nx Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(nxDefinition, "RetVal", "int32");
validate(nxDefinition);

%% C++ class method |ny| for C++ class |Image<short>| 
% C++ Signature: int Image<short>::ny()
nyDefinition = addMethod(Image_short_Definition, ...
    "int Image<short>::ny()", ...
    "MATLABName", "ny", ...
    "Description", "ny Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(nyDefinition, "RetVal", "int32");
validate(nyDefinition);

%% C++ class method |nz| for C++ class |Image<short>| 
% C++ Signature: int Image<short>::nz()
nzDefinition = addMethod(Image_short_Definition, ...
    "int Image<short>::nz()", ...
    "MATLABName", "nz", ...
    "Description", "nz Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(nzDefinition, "RetVal", "int32");
validate(nzDefinition);

%% C++ class method |pos| for C++ class |Image<short>| 
% C++ Signature: Vec3<double> Image<short>::pos()
posDefinition = addMethod(Image_short_Definition, ...
    "Vec3<double> Image<short>::pos()", ...
    "MATLABName", "pos", ...
    "Description", "pos Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(posDefinition, "RetVal", "clib.Image.Vec3_double_");
validate(posDefinition);

%% C++ class method |nFrames| for C++ class |Image<short>| 
% C++ Signature: int Image<short>::nFrames()
nFramesDefinition = addMethod(Image_short_Definition, ...
    "int Image<short>::nFrames()", ...
    "MATLABName", "nFrames", ...
    "Description", "nFrames Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(nFramesDefinition, "RetVal", "int32");
validate(nFramesDefinition);

%% C++ class method |nRepeats| for C++ class |Image<short>| 
% C++ Signature: int Image<short>::nRepeats()
nRepeatsDefinition = addMethod(Image_short_Definition, ...
    "int Image<short>::nRepeats()", ...
    "MATLABName", "nRepeats", ...
    "Description", "nRepeats Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(nRepeatsDefinition, "RetVal", "int32");
validate(nRepeatsDefinition);

%% C++ class method |nBytes| for C++ class |Image<short>| 
% C++ Signature: int Image<short>::nBytes()
nBytesDefinition = addMethod(Image_short_Definition, ...
    "int Image<short>::nBytes()", ...
    "MATLABName", "nBytes", ...
    "Description", "nBytes Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(nBytesDefinition, "RetVal", "int32");
validate(nBytesDefinition);

%% C++ class method |get| for C++ class |Image<short>| 
% C++ Signature: int16_t const * Image<short>::get()
getDefinition = addMethod(Image_short_Definition, ...
   "int16_t const * Image<short>::get()", ...
   "MATLABName", "get", ...
   "Description", "get Method of C++ class Image<short>."); % Modify help description values as needed.
defineOutput(getDefinition, "RetVal", "int16", ["nx","ny","nRepeats","nz","nFrames"]); % <MLTYPE> can be "int16", or "clib.array.Image.Short"
validate(getDefinition);

%% C++ class |SImage| with MATLAB name |clib.Image.SImage| 
SImageDefinition = addClass(libDef, "SImage", "MATLABName", "clib.Image.SImage", ...
    "Description", "clib.Image.SImage    Representation of C++ class SImage."); % Modify help description values as needed.

%% C++ class constructor for C++ class |SImage| 
% C++ Signature: SImage::SImage(char const * filename)
SImageConstructor1Definition = addConstructor(SImageDefinition, ...
   "SImage::SImage(char const * filename)", ...
   "Description", "clib.Image.SImage Constructor of C++ class SImage."); % Modify help description values as needed.
defineArgument(SImageConstructor1Definition, "filename", "string", "input", "nullTerminated"); % <MLTYPE> can be "clib.array.Image.Char","int8","string", or "char"
validate(SImageConstructor1Definition);

%% C++ class constructor for C++ class |SImage| 
% C++ Signature: SImage::SImage(SImage const & input1)
SImageConstructor2Definition = addConstructor(SImageDefinition, ...
    "SImage::SImage(SImage const & input1)", ...
    "Description", "clib.Image.SImage Constructor of C++ class SImage."); % Modify help description values as needed.
defineArgument(SImageConstructor2Definition, "input1", "clib.Image.SImage", "input");
validate(SImageConstructor2Definition);

%% Validate the library definition
validate(libDef);

end
