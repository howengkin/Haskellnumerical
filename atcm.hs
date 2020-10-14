-- This file contains all the Haskell programs for ATCM2017 paper on 
-- "Appreciating functions via functional programming"

add ::  Double -> Double -> Double
add x y = x + y

sub :: Double -> Double -> Double
sub x y = x - y

mult :: Double -> Double -> Double
mult x y = x*y

divd :: Double -> Double -> Double
divd x y = x/y

powr :: Double -> Double -> Double
powr x y = x**y

-- Factorial

fact :: Integer -> Integer
fact 0 = 1
fact n = n*fact (n-1)

-- Root finding
-- Bisection scheme

type Interval = (Double,Double)
type Function = (Double -> Double)

good :: Function -> Interval -> Bool
good f (a,b) = ((f a)*(f b) <= 0)

goodhalf :: Function -> Interval -> Interval
goodhalf f (a,b) = let c = (a+b)/2 in
                   if (good f (a,c)) then (a,c) 
				                     else (c,b)

bisectseq :: Function -> Interval -> [Interval]
bisectseq f (a,b) = (a,b): bisectseq f (goodhalf f (a,b))

bisect :: Double ->  Function -> Interval -> Double
bisect eps f (a,b) = relerrint eps (bisectseq f (a,b))

-- Relative error checker

relerrint :: Double -> [Interval] -> Double
relerrint eps ((a,b):xs)
 | abs(a/b - 1) < eps = a
 | otherwise          = relerrint eps xs 

relerr :: Double -> [Double] -> Double
relerr eps (a:b:xs)
 | abs(a/b - 1) < eps = a
 | otherwise          = relerr eps (b:xs)

-- Absolute error checker

abserr :: Double -> [Double] -> Double
abserr eps (a:b:xs)
 | abs(a-b) < eps = a
 | otherwise      = abserr eps (b:xs)
 
-- Linear interpolation scheme

linecut :: Function -> Interval ->  Double
linecut f (a,b) = (a*f(b)-b*f(a))/(f(b)-f(a))

linehalf :: Function -> Interval -> Interval
linehalf f (a,b) = let c = linecut f (a,b) in 
                  if (good f (a,c)) then (a,c) 
				                    else (c,b)

linehalfseq :: Function -> Interval -> [Interval]
linehalfseq f (a,b) = (a,b) : linehalfseq f (linehalf f (a,b))
									
linint :: Double -> Function -> Interval -> Double
linint eps f (a,b) = relerrint eps (linehalfseq f (a,b))

-- Fixed point iteration scheme
-- x(n+1) = f(x(n)), x(0) = a 

iterates :: Double -> Function -> [Double]
iterates a f = a : (iterates (f a) f)

fpi :: Double -> Double -> Function -> Double
fpi eps a f = relerr eps (iterates a f)

-- Example

f1 :: Function
f1 x = 1 + 1/x


-- Newton-Raphson iteration scheme
 
newtonit :: Function -> Function -> Function
newtonit f df x = x - (f x)/(df x)
 
newton :: Double -> Double -> Function -> Function -> Double
newton eps a f df = fpi eps a (newtonit f df)

-- Example
f2 :: Function
f2 x = x**3 + 2*(x**2) + 10*x - 20

df2 :: Function 
df2 x = 3*(x**2) + 4*x + 10

-- Numerical differentiation

secant :: Double -> Function -> Function 
secant h f x = (f (x+h) - f x)/h

secantseq :: Double -> Function -> Double -> [Double] 
secantseq h f x = (secant h f x) : secantseq (h/10) f x

easydiff :: Double -> Function -> Function
easydiff eps f x = relerr eps (secantseq 1 f x)

-- Example
f3 :: Function
f3 x = sin (x)

-- Euler's method for solving differential equation: y' = f(x,y)
-- step size h = (b-a)/n 
-- x_0 = a, x_n = b

type Function2 = Double -> Double -> Double

ynext :: Double -> Double -> Double -> Function2 -> Double
ynext h x y f = y + (h * (f x y))

eulerseq :: Double -> Double -> Double -> Function2 -> [Double]
eulerseq h x y f = w : (eulerseq h (x+h) w f)
    where w = ynext h x y f

euler :: Double -> Double -> Double -> Function2 -> Int -> Double
euler a b y0 f n = (eulerseq h a y0 f)!!n
                   where h = (b-a)/fromIntegral(n)

--memoise_auxeuler :: Int -> Function2 -> Double -> Double -> Double -> Double
--memoise_auxeuler = (map auxeuler [0..] !!)
--    where auxeuler 0 f h x y0 = y0
--          auxeuler m f h x y0 = (memoise_auxeuler (m-1) f h (x-h) y0) + h*(f (x-h) (memoise_auxeuler (m-1) f h (x-h) y0))
								
--auxeuler :: Function2 -> Double -> Double -> Int -> Double -> Double
--auxeuler f h x 0 y0 = y0
--auxeuler f h x m y0 = (auxeuler f h (x-h) (m-1) y0) + h*(f (x-h) (auxeuler f h (x-h) (m-1) y0))

--euler :: Function2 -> Double -> Double -> Int -> Double -> Double
--euler f a b n y0 = memoise_auxeuler n f ((b-a)/fromIntegral(n)) b y0  

-- Example
f4 :: Function2
f4 x y = 2*x + 2*y

-- Runge-Kutta 4 method

k1 :: Double -> Double -> Function2 -> Double
k1 x y f = f x y

k2 :: Double -> Double -> Double -> Function2 -> Double
k2 h x y f = f (0.5*h+x) (y+0.5*(k1 x y f)*h)

k3 :: Double -> Double -> Double -> Function2 -> Double
k3 h x y f = f (0.5*h+x) (y+0.5*(k2 h x y f)*h)

k4 :: Double -> Double -> Double -> Function2 -> Double
k4 h x y f = f (x+h) (y+(k3 h x y f) *h)

ksum :: Double -> Double -> Double -> Function2 -> Double
ksum h x y f = (k1 x y f) + 2*(k2 h x y f) + 2*(k3 h x y f) + (k4 h x y f)

rk4 :: Double -> Double -> Double -> Function2 -> Double
rk4 h x y f = y + (1/6)*h*(ksum h x y f)

rk4seq :: Double -> Double -> Double -> Function2 -> [Double]
rk4seq h x y f = w : rk4seq h (x+h) w f
    where w = rk4 h x y f

rk4meth :: Double -> Double -> Double -> Function2 -> Int -> Double
rk4meth a b y0 f n = (rk4seq h a y0 f)!!n 
                      where h = (b-a)/fromIntegral(n)

-- Example
f5 :: Function2
f5 x y = 3*exp(-x)-0.4*y

-- Numerical integration

trap :: Function -> Interval -> Double
trap f (a,b) = (f(a)+f(b))*(b-a)/2

zipadd :: [Double] -> [Double] -> [Double]
zipadd (x0:xs) (y0:ys) = (x0 + y0): zipadd xs ys

integralseq :: Function -> Interval -> [Double]
integralseq f (a,b) = (trap f (a,b)): zipadd (integralseq f (a,c)) (integralseq f (c,b))
                      where c = (a+b)/2

integrate :: Double -> Function -> Interval -> Double
integrate eps f (a,b) = relerr eps (drop 10 (integralseq f (a,b)))					  

-- Example

f6 :: Function
f6 x = x**2

fun1 :: Double -> Double
fun1 x = (exp(cos(2*pi*x)))*(cos(sin(2*pi*x)))

-- Complex numbers

type Complexr = (Double,Double)
type Complexp = (Double,Double)

fromptor :: Complexr -> Complexp
fromptor (r,t) = (r*cos(t),r*sin(t))

-- Real and imaginery parts

re :: Complexr -> Double
re (x,y) = x

im :: Complexr -> Double
im (x,y) = y

-- Complex arithmetic

addr :: Complexr -> Complexr -> Complexr
addr (x1,y1) (x2,y2) = (x1+x2,y1+y2)

minusr :: Complexr -> Complexr -> Complexr
minusr (x1,y1) (x2,y2) = (x1-x2,y1-y2)

multr :: Complexr -> Complexr -> Complexr
multr (x1,y1) (x2,y2) = (x1*x2-y1*y2,x1*y2+y1*x2)

divr :: Complexr -> Complexr -> Complexr
divr (x1,y1) (x2,y2) = ((x1*x2+y1*y2)/(x2**2+y2**2),
                        (y1*x2-x1*y2)/(x2**2+y2**2))

intpower :: Complexr -> Integer -> Complexr
intpower z 0 = (1,0)
intpower z n = multr z (intpower z (n-1))

conj :: Complexr -> Complexr
conj (x,y) = (x,-y)

-- Polynomial functions
polynom :: [Complexr] -> (Complexr -> Complexr)
polynom [] z = (0,0)
polynom (x:xs) z  = addr x (multr z (polynom xs z))  

-- Mobius transformation
mobius :: Complexr -> Complexr -> Complexr -> Complexr -> (Complexr -> Complexr)
mobius a b c d  z = divr (addr (multr a z) b) (addr (multr c z) d) 

-- complex cos(z), sin(z) and exp(z)

cosr :: Complexr -> Complexr
cosr (x,y) = ((cos x)*(cosh y),-(sin x)*(sinh y))

sinr :: Complexr -> Complexr
sinr (x,y) = ((sin x)*(cosh y),(cos x)*(sinh y))

expr :: Complexr -> Complexr
expr (x,y) = ((exp x)*(cos y),(exp x)*(sin y))

-- numerical differentiation using L-M method

build :: (Complexr -> Complexr) -> Integer -> Double -> Double -> (Double -> Double)
build f n a r = \t -> ((cos (2*pi*fromIntegral(n)*t))*re(f(a+r*cos(2*pi*t),r*sin(2*pi*t)))) 

nderiv :: (Complexr -> Complexr) -> Integer -> Double -> Double -> Double -> Double
nderiv f n a r eps = k*(integrate eps (build f n a r) (0,1)) 
					where k = (2*fromIntegral(fact n))/(r^n)
					
testfn :: Complexr -> Complexr
testfn z = divr (expr z) (addr (intpower (sinr z) 3) (intpower (cosr z) 3))