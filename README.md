This is a reference code for the poster `Guided Path Tracing Using Clustered Virtual Point Lights' in SIGGRAPH Asia 2015. 

The code is made to educate myself about physically based rendering a long time ago, with inspiration drawn from PBRT and Mitsuba. 
Please feel free to use it for your research. 

## Data

The code includes Breakfast, an example scene originally modelled by Wig42 on [BlendSwap](https://www.blendswap.com/blends/view/75431). 
The painting on the wall is digitized from an acryllic painting by my little sister [Kim](http://phuonghua.wixsite.com/phuong). 

## Compilation

To compile, please use Visual Studio 2013 Community edition. 
Unpack `precompiled.zip` which contains all required dependencies to the source code folder. It also contains a precompiled Windows x64 binary for the renderer.

To run, go to `dist` folder, and execute

	renzoku.exe scenes/breakfast/breakfast_metro.json

The images will be saved to the `data` folder in EXR format. 	

## Citation

If you find the code useful for your research, please cite the poster: 

	@inproceedings{hua-guided-sa15,
	 author = {Hua, Binh-Son and Low, Kok-Lim},
	 title = {Guided Path Tracing Using Clustered Virtual Point Lights},
	 booktitle = {SIGGRAPH Asia 2015 Posters},
	 year = {2015},
	 publisher = {ACM}
	} 

## Dependencies

For the convenience of compilation, this code is distributed with the following open source libraries:

* Boost 

* Eigen

* FreeGLUT

* FreeImage

* GLEW 

* LibImportance by Jiri Vorba under GNU General Public License v3. 

* libjson 

* Morton sort by NVIDIA

* OpenEXR

## License 

This code is released under MIT license. 

Copyright 2012-2018, Binh-Son Hua.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.