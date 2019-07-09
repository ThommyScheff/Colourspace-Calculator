# Colourspace-Calculator
A Colourspace Calculator which can convert every possible colourspace written in the code into eachother
Many thanks to every source used!
Sources included:

https://www.rapidtables.com
https://www.adobe.com/digitalimag/pdfs/AdobeRGB1998.pdf
http://www.brucelindbloom.com
http://www.easyrgb.com/en/math.php#text2

Every Class included:

Adobe RGB 1998 (RGB)
Standard RGB (sRGB) 
HSL in colourspace sRGB (Hue, Saturation, Lightness)
HSV in colourspace sRGB (Hue, Saturation, Value)
CIEXYZ
CIELab
CIELuv

Example input:

>>> col = RGB([1,0,0])
>>> XYZ = col.toCIEXYZ()
>>> XYZ.values
>>> [0.57667, 0.29734, 0.02703]
A successful conversion from RGB to CIEXYZ.
