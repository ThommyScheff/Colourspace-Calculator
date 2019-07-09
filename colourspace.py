#!/usr/bin/env python
"""
Colour Calculator which can convert every possible colourspace into eachother.
"""
__author__ = "Thommy Scheffner"
__copyright__ = "Copyright 2019, Goethe University for DocStrings and for 'do not modify this method' methods"
__credits__ = ["Thommy Scheffner"]

__license__ = "GPL"
__version__ = "3.7.3"
__maintainer__ = "Thommy Scheffner"
__email__ = "itzthommy@gmx.de"
__status__ = "Completed"



class colour():
    def __init__(self, values):
        self.values = values

    def __str__(self):
        output = ''
        for entry in self.values:
            output += '%.2f ' % entry
        return output


"""We always assume the reference white to be illuminant D65 with a 2° observer.
Always use the values *in the standard*, not the intended values.
"""

class RGB(colour):
    """A class for handling RGB values. Every colour lies in [0, 1].
       We are using the reference white D65 and the primary colors
       of the Adobe RGB Space from 1998. You may assume gamma to be 2.2
       Source for matrix values D65 and 2° which is used
       for every RGB -> XYZ, XYZ - RGB conversion is the Adobe RGB 1998 Paper
       """
    # Source: https://www.adobe.com/digitalimag/pdfs/AdobeRGB1998.pdf

    def __init__(self,values):
        self.values = values
        self.r = values[0]
        self.g = values[1]
        self.b = values[2]

    def tosRGB(self):
        """Converting RGB values to sRGB values with an intermediate step to
        the XYZ colourspace
        RGB -> XYZ -> sRGB"""
        newXYZ = RGB.toCIEXYZ(self)

        newsRGB = newXYZ.tosRGB()

        return newsRGB

    def toHSL(self):
        """Returns HSL values from given RGB values. Hue [0,359],
        Saturation [0,1] and Lightness [0,1]
        RGB -> XYZ -> sRGB -> HSL"""
        # Source for equations:
        # https://www.rapidtables.com/convert/color/rgb-to-hsl.html
        # HSL is in colourspace sRGB, so RGB values convert to sRGB values
        newsRGB = RGB.tosRGB(self)
        maximum = max(newsRGB.values)
        minimum = min(newsRGB.values)
        self.r = newsRGB.values[0]
        self.g = newsRGB.values[1]
        self.b = newsRGB.values[2]

        delta = maximum - minimum

        if maximum == minimum:
            h = 0
        elif maximum == self.r:
            h = ((self.g - self.b) / delta)
            h *= 60

        elif maximum == self.g:
            h = (2+ ((self.b - self.r) / delta))
            h *= 60

        elif maximum == self.b:
            h = (4 + ((self.r - self.g) / delta))
            h *= 60
        else:
            pass

        if h < 0:
            h += 360

        L = (maximum + minimum) / 2

        if maximum == 0:
            s = 0
        elif minimum == 1:
            s = 0
        else:
            betrag = 2 * L - 1

            if betrag < 0:
                betrag += betrag * -2

            s = delta / (1 - betrag)

        return HSL([h,s,L])

    def toHSV(self):
        """Returns HSV values from given RGB values. Hue, Saturation and
        Value"""
        # Source for equations:
        # https://www.rapidtables.com/convert/color/rgb-to-hsv.html
        # HSV is in colourspace sRGB, so RGB values convert to sRGB values
        newsRGB = RGB.tosRGB(self)
        maximum = max(newsRGB.values)
        minimum = min(newsRGB.values)
        self.r = newsRGB.values[0]
        self.g = newsRGB.values[1]
        self.b = newsRGB.values[2]
        
        delta = maximum - minimum

        if maximum == minimum:
            h = 0
        elif maximum == self.r:
            h = ((self.g - self.b) / delta)
            h *= 60

        elif maximum == self.g:
            h = (2+ ((self.b - self.r) / delta))
            h *= 60

        elif maximum == self.b:
            h = (4 + ((self.r - self.g) / delta))
            h *= 60
        else:
            pass

        if h < 0:
            h += 360
            
        if maximum == 0:
            s = 0
        else:
            s = delta / maximum

        v = maximum

        return HSV([h,s,v])

    def toCIEXYZ(self):
        """Returns XYZ values from given RGB values. Reference white D65, 2°
        and primary colours for RGB from the Adobe RGB Paper 1998
        Gamma is 2.2"""
        # Source: https://www.adobe.com/digitalimag/pdfs/AdobeRGB1998.pdf
        
        r = self.r ** 2.2
        g = self.g ** 2.2
        b = self.b ** 2.2
        
        # Adobe RGB 1998 Paper
        x = r * 0.57667 + g * 0.18556 + b * 0.18823
        y = r * 0.29734 + g * 0.62736 + b * 0.07529
        z = r * 0.02703 + g * 0.07069 + b * 0.99134

        return CIEXYZ([x,y,z])

    def toCIELab(self):
        """Returns Lab values from given RGB values. RGB -> XYZ -> Lab.
        Reference white D65, 2° and primary colours from Adobe RGB 1998
        (Bruce Lindbloom). Gamma is 2.2"""
        # Convert RGB to XYZ so we can use those values for CIELab

        newXYZ = RGB.toCIEXYZ(self)
        x = newXYZ.values[0]
        y = newXYZ.values[1]
        z = newXYZ.values[2]
        
        newLab = newXYZ.toCIELab()

        return newLab

    def toCIELuv(self):
        """Returns Luv values from given RGB values. RGB -> XYZ -> Luv.
        Reference white D65, 2° and primary colours from Adobe RGB 1998
        (Bruce Lindbloom). Gamma is 2.2"""
        # Convert RGB to XYZ so we can use those values for CIELuv
        newXYZ = RGB.toCIEXYZ(self)
        x = newXYZ.values[0]
        y = newXYZ.values[1]
        z = newXYZ.values[2]

        newLuv = newXYZ.toCIELuv()

        return newLuv

class sRGB(colour):
    """A class for handling RGB values. Every colour lies in [0, 1].
       We are using the reference white D65 and the primary colors
       of the sRGB Color Space. Gamma companding follows the sRGB standard"""

    def __init__(self, values):
        self.values = values
        self.r = values[0]
        self.g = values[1]
        self.b = values[2]

    def toRGB(self):
        """Returns RGB values from given sRGB values. We are using the
        reference white D65 and the primary colors of the RGB Color Space.
        Gamma is 2.2 and R,G,B : [0,1].
        sRGB -> CIEXYZ -> Adobe RGB"""
        
        newXYZ = sRGB.toCIEXYZ(self)

        newRGB = newXYZ.toRGB()

        return newRGB
    
    def toHSL(self):
        """Returns HSL values from given sRGB values. Hue [0,359],
        Saturation [0,1] and Lightness [0,1]
        sRGB -> HSL"""
        # Source for equations:
        # https://www.rapidtables.com/convert/color/rgb-to-hsl.html
        maximum = max(self.values)
        minimum = min(self.values)
        delta = maximum - minimum

        if maximum == minimum:
            h = 0
        elif maximum == self.r:
            h = ((self.g - self.b) / delta)
            h *= 60

        elif maximum == self.g:
            h = (2+ ((self.b - self.r) / delta))
            h *= 60

        elif maximum == self.b:
            h = (4 + ((self.r - self.g) / delta))
            h *= 60
        else:
            pass

        if h < 0:
            h += 360

        L = (maximum + minimum) / 2

        if maximum == 0:
            s = 0
        elif minimum == 1:
            s = 0
        else:
            betrag = 2 * L - 1

            if betrag < 0:
                betrag += betrag * -2

            s = delta / (1 - betrag)

        return HSL([h,s,L])
    
    def toHSV(self):
        """Returns HSV values from given sRGB values. Hue, Saturation and
        Value"""
        # Source for equations:
        # https://www.rapidtables.com/convert/color/rgb-to-hsv.html
        maximum = max(self.values)
        minimum = min(self.values)
        delta = maximum - minimum

        if maximum == minimum:
            h = 0
        elif maximum == self.r:
            h = ((self.g - self.b) / delta)
            h *= 60

        elif maximum == self.g:
            h = (2+ ((self.b - self.r) / delta))
            h *= 60

        elif maximum == self.b:
            h = (4 + ((self.r - self.g) / delta))
            h *= 60
        else:
            pass

        if h < 0:
            h += 360
            
        if maximum == 0:
            s = 0
        else:
            s = delta / maximum

        v = maximum

        return HSV([h,s,v])
    
    def toCIEXYZ(self):
        """Returns XYZ values from given sRGB values. Reference white D65, 2°
        and primary colours for sRGB from: 
        http://www.brucelindbloom.com/
        Gamma is 2.4"""
        if self.r <= 0.04045:
            r = self.r / 12.92
        else:
            r = ((self.r+0.055)/1.055)**2.4
        if self.g <= 0.04045:
            g = self.g / 12.92
        else:
            g = ((self.g+0.055)/1.055)**2.4
        if self.b <= 0.04045:
            b = self.b / 12.92
        else:
            b = ((self.b+0.055)/1.055)**2.4
            
        # -> Some Common RGB Working Space Matrices -> sRGB to XYZ[M]
        x = r * 0.4124564 + g * 0.3575761 + b * 0.1804375
        y = r *  0.2126729 + g * 0.7151522 + b * 0.0721750
        z = r * 0.0193339 + g * 0.1191920 + b * 0.9503041

        return CIEXYZ([x,y,z])

    def toCIELab(self):
        """Returns Lab values from given sRGB values. sRGB -> XYZ -> Lab.
        Reference white D65, 2° and primary colours from Adobe RGB 1998.
        Gamma is 2.2"""
        newXYZ = sRGB.toCIEXYZ(self)
        x = newXYZ.values[0]
        y = newXYZ.values[1]
        z = newXYZ.values[2]

        newLab = newXYZ.toCIELab()
        
        return newLab
    
    def toCIELuv(self):
        """Returns Luv values from given RGB values. RGB -> XYZ -> Luv.
        Reference white D65, 2° and primary colours from Adobe RGB 1998
        (Bruce Lindbloom). Gamma is 2.2"""
        # Convert sRGB to XYZ so we can use those values for CIELuv
        newXYZ = sRGB.toCIEXYZ(self)
        x = newXYZ.values[0]
        y = newXYZ.values[1]
        z = newXYZ.values[2]

        newLuv = newXYZ.toCIELuv()

        return newLuv

class HSL(colour):
    """HSL colour model for the sRGB colour space.
    Use the same assumptions as when working with sRGB.
    Greyscale colours have a H of 0, rather than undefined.
    """

    def __init__(self, values):
        self.values = values
        self.H = values[0]
        self.S = values[1]
        self.L = values[2]

    def toRGB(self):
        """Returns RGB values from given HSL values. We are using the
        reference white D65 and the primary colors of the RGB Color Space.
        Gamma is 2.2 and R,G,B : [0,1].
        HSL -> sRGB -> CIEXYZ -> Adobe RGB"""
        #Source: https://www.rapidtables.com/convert/color/hsl-to-rgb.html

        newsRGB = HSL.tosRGB(self)
        
        newXYZ = newsRGB.toCIEXYZ()
        
        newRGB = newXYZ.toRGB()

        return newRGB

    def tosRGB(self):
        """Converting HSL values to sRGB values.
        Source:
        https://www.rapidtables.com/convert/color/hsl-to-rgb.html
        HSL -> sRGB"""
        
        # betrag = absolute value of an equation -> no negative result
        betrag = (2 * self.L - 1)
        if betrag < 0:
            betrag += betrag * -2
    
        betrag2 = (self.H / 60) % 2 - 1
        if betrag2 < 0:
            betrag2 += betrag2 * -2
            
        C = (1 - betrag) * self.S
        X = C * (1 - betrag2)
        m = self.L - C / 2

        if self.H >= 0 and self.H < 60:
            R = C
            G = X
            B = 0
        elif self.H >= 60 and self.H < 120:
            R = X
            G = C
            B = 0
        elif self.H >= 120 and self.H < 180:
            R = 0
            G = C
            B = X
        elif self.H >= 180 and self.H < 240:
            R = 0
            G = X
            B = C
        elif self.H >= 240 and self.H < 300:
            R = X
            G = 0
            B = C
        elif self.H >= 300 and self.H < 360:
            R = C
            G = 0
            B = X
        else:
            pass

        R = R + m
        G = G + m
        B = B + m
        
        return sRGB([R,G,B])

    def toHSV(self):
        """Returns HSV values from given RGB values. Hue, Saturation and
        Value"""
        # Source for equations:
        # https://www.rapidtables.com/convert/color/hsl-to-rgb.html
        # Conversion HSL -> sRGB -> HSV

        # HSL -> sRGB now:
        newsRGB = HSL.tosRGB(self)
        # HSL -> sRGB converted:
        
        # sRGB -> HSV now:
        newHSV = newsRGB.toHSV()
        # sRGB -> HSV converted:
        
        # Source for equations:
        # https://www.rapidtables.com/convert/color/rgb-to-hsv.html

        return newHSV

    def toCIEXYZ(self):
        """Returns XYZ values from given HSL values. Reference white D65, 2°
        and primary colours for sRGB from: 
        http://www.brucelindbloom.com/
        Gamma is 2.4
        Converting HSL -> sRGB -> CIEXYZ"""
        # https://www.rapidtables.com/convert/color/hsl-to-rgb.html
        # betrag = absolute value of an equation -> no negative result
        # Conversion HSL -> sRGB -> XYZ

        newsRGB = HSL.tosRGB(self)
        # HSL -> sRGB converted

        # sRGB -> CIEXYZ now
        newXYZ = newsRGB.toCIEXYZ()

        return newXYZ

    def toCIELab(self):
        """Returns Lab values from given HSL values.
        HSL -> sRGB -> CIEXYZ -> CIELab.
        Reference white D65, 2° and primary colours from Adobe RGB 1998.
        Gamma is 2.2"""
        # https://www.rapidtables.com/convert/color/hsl-to-rgb.html
        # betrag = absolute value of an equation -> no negative result
        # Conversion HSL -> sRGB -> XYZ -> Lab

        newsRGB = HSL.tosRGB(self)

        newXYZ = newsRGB.toCIEXYZ()

        newLab = newXYZ.toCIELab()

        return newLab

    def toCIELuv(self):
        # https://www.rapidtables.com/convert/color/hsl-to-rgb.html
        # betrag = absolute value of an equation -> no negative result
        # Conversion HSL -> sRGB -> XYZ -> Luv

        newsRGB = HSL.tosRGB(self)

        newXYZ = newsRGB.toCIEXYZ()

        newLuv = newXYZ.toCIELuv()

        return newLuv

class HSV(colour):
    """HSV colour model for the sRGB colour space.
    Use the same assumptions as when working with sRGB.
    Greyscale colours have a H of 0, rather than undefined.
    """

    def __init__(self, values):
        self.values = values
        self.H = values[0]
        self.S = values[1]
        self.V = values[2]
        
    def toRGB(self):
        """Returns RGB values from given HSL values. We are using the
        reference white D65 and the primary colors of the RGB Color Space.
        Gamma is 2.2 and R,G,B : [0,1].
        HSV -> sRGB -> CIEXYZ -> Adobe RGB"""

        newsRGB = HSV.tosRGB(self)

        newXYZ = newsRGB.toCIEXYZ()

        newRGB = newXYZ.toRGB()
        
        return newRGB

    def tosRGB(self):
        """Converting HSV values to sRGB values.
        Source:
        https://www.rapidtables.com/convert/color/hsv-to-rgb.html
        HSV -> sRGB"""
        #betrag = absolute value of equation, no negative outcome
        
        betrag = (self.H / 60) % 2 - 1
        if betrag < 0:
            betrag += betrag * -2

        C = self.V * self.S
        X = C * (1 - betrag)
        m = self.V - C

        if self.H >= 0 and self.H < 60:
            R = C
            G = X
            B = 0
        elif self.H >= 60 and self.H < 120:
            R = X
            G = C
            B = 0
        elif self.H >= 120 and self.H < 180:
            R = 0
            G = C
            B = X
        elif self.H >= 180 and self.H < 240:
            R = 0
            G = X
            B = C
        elif self.H >= 240 and self.H < 300:
            R = X
            G = 0
            B = C
        elif self.H >= 300 and self.H < 360:
            R = C
            G = 0
            B = X
        else:
            pass

        R = R + m
        G = G + m
        B = B + m

        return sRGB([R,G,B])

    def toHSL(self):
        """Returns HSL values from given HSV values. Hue [0,359],
        Saturation [0,1] and Lightness [0,1]
        HSV -> HSL"""
        L = (2 - self.S) * self.V / 2
        if L != 0:
            if L == 1:
                s = 0
            elif L < 0.5:
                S = self.S * self.V / (L * 2)
            else:
                S = self.S * self.V / (2 - L * 2)

        return colour([self.H, S, L])


    def toCIEXYZ(self):
        """Returns XYZ values from given HSV values. Reference white D65, 2°
        and primary colours for sRGB from: 
        http://www.brucelindbloom.com/
        Gamma is 2.4
        Converting HSV -> sRGB -> CIEXYZ"""
        # HSV -> sRGB now:
        newsRGB = HSV.tosRGB(self)
        # HSV -> sRGB completed
        
        # sRGB -> XYZ now:
        newXYZ = sRGB([r,g,b]).toCIEXYZ()
        # sRGB -> XYZ completed
        
        return newXYZ

    def toCIELab(self):
        """Returns Lab values from given HSL values.
        HSV -> sRGB -> CIEXYZ -> CIELab.
        Reference white D65, 2° and primary colours from Adobe RGB 1998.
        Gamma is 2.2"""
        # HSV - sRGB now:
        newsRGB = HSV.tosRGB(self)
        r = newsRGB.values[0]
        g = newsRGB.values[1]
        b = newsRGB.values[2]
        # HSV -> sRGB completed
        
        # sRGB -> XYZ now:
        newXYZ = sRGB([r,g,b]).toCIEXYZ()
        x = newXYZ.values[0]
        y = newXYZ.values[1]
        z = newXYZ.values[2]
        # sRGB -> XYZ completed
        
        # XYZ -> Lab now:
        newLab = CIEXYZ([x,y,z]).toCIELab()

        return newLab

    def toCIELuv(self):
        # HSV -> sRGB now:
        newsRGB = HSV.tosRGB(self)
        r = newsRGB.values[0]
        g = newsRGB.values[1]
        b = newsRGB.values[2]
        # HSV -> sRGB completed
        
        # sRGB -> XYZ now:
        newXYZ = sRGB([r,g,b]).toCIEXYZ()
        x = newXYZ.values[0]
        y = newXYZ.values[1]
        z = newXYZ.values[2]
        # sRGB -> XYZ completed

        # XYZ -> Luv now:
        newLuv = CIEXYZ([x,y,z]).toCIELuv()

        return newLuv
        

class CIEXYZ(colour):
    """CIE XYZ colour space. Every value lies within [0, 1]"""

    def __init__(self, values):
        self.values = values
        self.x = values[0]
        self.y = values[1]
        self.z = values[2]

    def toRGB(self):
        """Returns RGB values from given XYZ values. We are using the
        reference white D65 and the primary colors of the RGB Color Space.
        Gamma is 2.2 and R,G,B : [0,1].
        HSL -> sRGB -> CIEXYZ -> Adobe RGB"""
        # Matrix according to source:
        # Adobe RGB Paper 1998
        r = self.x * 2.04159 + self.y * -0.56501 + self.z * -0.34473
        g = self.x * -0.96924 + self.y * 1.87597 + self.z * 0.04156
        b = self.x * 0.01344 + self.y * -0.11836 + self.z * 1.01517

        #gamma = 2.2 for Adobe RGB with D65 and 2
        R = r**(1/2.2)
        G = g**(1/2.2)
        B = b**(1/2.2)

        return RGB([R,G,B])

    def tosRGB(self):
        """Converting XYZ values to sRGB values.
        Source:
        http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html 
        CIEXYZ -> sRGB"""
        # Matrix according to source:
        # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html        
        r = self.x * 3.2404542 + self.y * -1.5371385 + self.z * -0.4985314
        g = self.x * -0.9692660 + self.y * 1.8760108 + self.z * 0.0415560
        b = self.x * 0.0556434 + self.y * -0.2040259 + self.z * 1.0572252

        # Gamma = 2.4 for Standard RGB
        R = r**(1/2.4)
        G = g**(1/2.4)
        B = b**(1/2.4)

        # sRGB Companding according to source:
        # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
        if r <= 0.0031308:
            R = 12.92 * r
        else:
            R = (1.055*(r**(1/2.4))) -0.055
        if g <= 0.0031308:
            G = 12.92 * g
        else:
            G = (1.055*(g**(1/2.4))) -0.055
        if b <= 0.0031308:
            B = 12.92 * b
        else:
            B = (1.055*(b**(1/2.4))) -0.055

        return sRGB([R,G,B])

    def toHSL(self):
        """Returns HSL values from given sRGB values. Hue [0,359],
        Saturation [0,1] and Lightness [0,1]
        XYZ -> sRGB -> HSL"""
        ###XYZ conversion to sRGB because HSL is in colourspace sRGB
        ##matrix according to source:
        #http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html        

        # XYZ -> sRGB now:
        newsRGB = CIEXYZ.tosRGB(self)
        R = newsRGB.values[0]
        G = newsRGB.values[1]
        B = newsRGB.values[2]
        # XYZ -> sRGB completed

        # sRGB -> HSL now:
        newHSL = newsRGB.toHSL()
        # sRGB -> HSL completed

        return newHSL

    def toHSV(self):
        """Returns HSV values from given XYZ values. Hue [0,359],
        Saturation [0,1] and Value [0,1]"""
        #XYZ conversion to sRGB because HSV is in colourspace sRGB
        #matrix according to source:
        #http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
        
        # XYZ -> sRGB now:
        newsRGB = CIEXYZ.tosRGB(self)
        # XYZ -> sRGB completed

        # sRGB -> HSV now:
        newHSV = newsRGB.toHSV()
        # sRGB -> HSV completed

        return newHSV

    def toCIELab(self):
        """Returns Lab values from given XYZ values.
        CIEXYZ -> CIELab.
        Reference white D65, 2° and primary colours from Adobe RGB 1998.
        Gamma is 2.2"""
        # Source for conversion:
        # http://www.brucelindbloom.com XYZ -> Lab 
        k = 903.3
        e = 0.008856
  
        xw = 0.95047
        yw = 1.00
        zw = 1.08883

        var_X = self.x / xw
        var_Y = self.y / yw
        var_Z = self.z / zw

        if  (var_X > e):
            fx = pow(var_X, (1/3))
        else:
            fx = (k* var_X + 16) / 116
            
        if  (var_Y > e):
            fy = pow(var_Y, (1/3))
            #var_Y = var_Y **(1/3)
        else:
            fy = (k* var_Y + 16) / 116
            
        if  (var_Z > e):
            fz = pow(var_Z, (1/3))
            #var_Z = var_Z **(1/3)
        else:
            fz = (k* var_Z + 16) / 116

        L = 116 * fy - 16
        a = 500 * (fx - fy)
        b = 200 * (fy - fz)

        return CIELab([L,a,b])
    
    def toCIELuv(self):
        #reference white values, source: Adobe RGB 1998
        #
        x = self.x
        y = self.y
        z = self.z

        xw = 0.95047
        yw = 1
        zw = 1.08883

        yr = y / yw

        utest = (4*x)/(x+15*y+3*z)
        vtest = (9*y)/(x+15*y+3*z)
        ur = (4*xw)/(xw+15*yw+3*zw)
        vr = (9*yw)/(xw+15*yw+3*zw)


        if yr > 0.008856:
            L = (116*(yr**(1/3)))-16
        else:
            L = 903.3 * yr

        u = 13*L*(utest - ur)
        v = 13*L*(vtest - vr)

        return colour([L,u,v])

class CIELab(colour):
    """CIELab colour space"""

    def __init__(self, values):
        self.values = values
        self.L = values[0]
        self.a = values[1]
        self.b = values[2]

    def toRGB(self):
        """Returns RGB values from given Lab values. We are using the
        reference white D65 and the primary colors of the RGB Color Space.
        Gamma is 2.2 and R,G,B : [0,1].
        Lab -> CIEXYZ -> Adobe RGB"""
        
        # Lab -> XYZ now:
        newXYZ = CIELab.toCIEXYZ(self)
        # Lab -> XYZ completed

        # XYZ -> RGB now:
        newRGB = newXYZ.toRGB()
        # XYZ -> RGB completed

        return newRGB
        
    def tosRGB(self):
        """Converting Lab values to sRGB values.
        Source:
        http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html 
        CIELab -> CIEXYZ -> sRGB"""

        # Lab -> XYZ now:
        newXYZ = CIELab.toCIEXYZ(self)
        # Lab -> XYZ completed

        # XYZ -> sRGB now:
        newsRGB = newXYZ.tosRGB()
        # XYZ -> sRGB completed

        return newsRGB
    
    def toHSL(self):
        """Returns HSL values from given sRGB values. Hue [0,359],
        Saturation [0,1] and Lightness [0,1]
        Lab -> XYZ -> sRGB -> HSL"""

        # Lab -> XYZ now:
        newXYZ = CIELab.toCIEXYZ(self)
        # Lab -> XYZ completed

        # XYZ -> sRGB now:
        newsRGB = newXYZ.tosRGB()
        # XYZ -> sRGB completed

        # sRGB -> HSL now:
        newHSL = newsRGB.toHSL()
        # sRGB -> HSL completed

        return newHSL

    def toHSV(self):
        """Returns HSV values from given Lab values. Hue [0,359],
        Saturation [0,1] and Value [0,1]
        Lab -> XYZ -> sRGB -> HSV"""
        # Lab -> XYZ now:
        newXYZ = CIELab.toCIEXYZ(self)
        # Lab -> XYZ completed

        # XYZ -> sRGB now:
        newsRGB = newXYZ.tosRGB()
        # XYZ -> sRGB completed

        # sRGB -> HSV now:
        newHSV = newsRGB.toHSV()
        # sRGB -> HSV completed

        return newHSV
    
    def toCIEXYZ(self):
        """Returns XYZ values from given Lab values. Reference white D65, 2°
        and primary colours for sRGB from: 
        http://www.brucelindbloom.com/
        Gamma is 2.4
        Converting CIELab -> CIEXYZ"""
        # Source for equation:
        # http://www.brucelindbloom.com/ -> Lab to XYZ
        k = 903.3
        e = 0.008856
        xw = 0.95047
        yw = 1
        zw = 1.08883

        fy = (self.L + 16) / 116
        fx = self.a / 500 + fy
        fz = fy - self.b / 200

        if fx**3 > e:
            xr = fx**3
        else:
            xr = (116*fx-16) / k
            
        if self.L > (e * k):
            yr = ((self.L + 16) / 116)**3
        else:
            yr = self.L / k
            
        if fz**3 > e:
            zr = fz**3
        else:
            zr = (116*fz-16) / k

        X = xr * xw
        Y = yr * yw
        Z = zr * zw

        return CIEXYZ([X,Y,Z])
        
    def toCIELuv(self):
        
        # Lab -> XYZ now:
        newXYZ = CIELab.toCIEXYZ(self)
        # Lab -> XYZ completed

        # XYZ -> Luv now:
        newLuv = newXYZ.toCIELuv()
        # XYZ -> Luv now:
        
        return newLuv


class CIELuv(colour):
    """CIELuv colour space"""
    

    def __init__(self, values):
        self.values = values
        self.L = values[0]
        self.u = values[1]
        self.v = values[2]

    def toRGB(self):
        """Returns RGB values from given Luv values. We are using the
        reference white D65 and the primary colors of the RGB Color Space.
        Gamma is 2.2 and R,G,B : [0,1].
        Luv -> CIEXYZ -> Adobe RGB"""
        
        # Luv -> XYZ now:
        newXYZ = CIELuv.toCIEXYZ(self)
        # Luv -> XYZ completed

        # XYZ -> RGB now:
        newRGB = newXYZ.toRGB()
        # XYZ -> RGB completed

        return newRGB
        
    def tosRGB(self):
        """Converting Luv values to sRGB values.
        Source:
        http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html 
        CIELuv -> CIEXYZ -> sRGB"""

        # Luv -> XYZ now:
        newXYZ = CIELuv.toCIEXYZ(self)
        # Luv -> XYZ completed

        # XYZ -> sRGB now:
        newsRGB = newXYZ.tosRGB()
        # XYZ -> sRGB completed

        return newsRGB
    
    def toHSL(self):
        """Returns HSL values from given sRGB values. Hue [0,359],
        Saturation [0,1] and Lightness [0,1]
        Luv -> XYZ -> sRGB -> HSL"""

        # Luv -> XYZ now:
        newXYZ = CIELuv.toCIEXYZ(self)
        # Luv -> XYZ completed

        # XYZ -> sRGB now:
        newsRGB = newXYZ.tosRGB()
        # XYZ -> sRGB completed

        # sRGB -> HSL now:
        newHSL = newsRGB.toHSL()
        # sRGB -> HSL completed

        return newHSL

    def toHSV(self):
        """Returns HSV values from given Luv values. Hue [0,359],
        Saturation [0,1] and Value [0,1]
        Luv -> XYZ -> XYZ -> sRGB -> HSV"""
        # Luv -> XYZ now:
        newXYZ = CIELuv.toCIEXYZ(self)
        # Luv -> XYZ completed

        # XYZ -> sRGB now:
        newsRGB = newXYZ.tosRGB()
        # XYZ -> sRGB completed

        # sRGB -> HSV now:
        newHSV = newsRGB.toHSV()
        # sRGB -> HSV completed

        return newHSV
    
    def toCIEXYZ(self):
        """Returns XYZ values from given Luv values. Reference white D65, 2°
        and primary colours for sRGB from: 
        http://www.brucelindbloom.com/
        Gamma is 2.4
        Converting CIELuv -> CIEXYZ"""
        # Source for Code:
        # http://www.easyrgb.com/en/math.php#text2
        # CIELuv -> XYZ
        k = 903.3
        e = 0.008856
        xw = 0.95047
        yw = 1
        zw = 1.08883

        y = (self.L + 16) / 116
        if (y**3) > e :
            y = y**3
        else:
            y= (y - 16 / 116) / 7.787

        u0 = (4 * xw) / (xw + (15 * yw) + (3 * zw))
        v0 = (9 * yw) / (xw + (15 * yw) + (3 * zw))

        u = self.u / ( 13 * self.L ) + u0
        v = self.v / ( 13 * self.L ) + v0

        Y = y * 100
        X =  -(9 * Y * u) / ((u - 4) * v - u * v)
        Z = (9 * Y -(15 * v * Y) - (v * X)) / (3 * v)

        X /= 100
        Y /= 100
        Z /= 100
        
        return CIEXYZ([X,Y,Z])
    
    def toCIELab(self):
        """Returns Lab values from given Luv values.
        CIELuv -> CIEXYZ -> CIELab.
        Reference white D65, 2° and primary colours from Adobe RGB 1998.
        Gamma is 2.2"""
        
        # Luv -> XYZ now:
        newXYZ = CIELuv.toCIEXYZ(self)
        # Luv -> XYZ completed

        newLab = newXYZ.toCIELab()

        return newLab


