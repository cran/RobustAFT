.dFvGet <- function()
{
    if(!exists(".dFv", envir=.GlobalEnv))
        dfvals()
    get(".dFv", envir=.GlobalEnv)
}

.dFvSet <- function(def)
{
    assign(".dFv", def, envir=.GlobalEnv, inherits=TRUE)
}

"dfvals" <-
function() {
assign(".dFv",list(
#
# Default values for LCMAIN / RGMAIN
#
tlo = 1.e-3 , gma = 1.0    , mxs = 1 , mxt = 30 , ntm = 0,
tua = 1.e-6 , tlu = 1.e-11 , iop = 2 , ix1 = 1  , iy1 = 1,
ic1 = 1     , ini = 2      , isr = 1 , itc = 1  , icn = 1,
alf = 0.025 , ccc = 1.345  , upr = 10.0,tli = 1.e-4,isq = 1,
isg = 1     , ite = 1      , itw = 0,
#
# Default values for WGMAIN
#
mxf = 150  , mxn = 50 , mxg = 50 , iwg = 2,
apr = 1.25 , icv = 1  , xfd = 1.0,
#
# Default values for KWMAIN
#
ia1 = 1  , ia2 = 1  , fff = 0.0 , ff1 = 0.0 , fu1 = 1.0 , fb1 = 1.0,
#
# Default values for AEMAIN
#
ial = 1  , mxe = 150,
#
# Default values for TSMAIN
#
fct = 0.0 , ffc  = 0.0,
#
# Default values for HBMAIN
#
ial = 1  , mxe = 150,
tls = 1.e-3  , tlr = 1.e-3  , msx = 30 , ik1 = 2 , ipt = 1,
isd = 1313   , ich = 1,
#
# Default values for CVMAIN
#
esp = 0.1  , ilc = 1 , aa2 = 0 , bb2  = 9, em=1.345, cr=2,  enu=1,
#
# Default values for MXMAIN
#
tlv = -6.9078  , tlm = 1.e-3 , ith = 1 , ilm = 1 , ddd = 1.345 ,
#
# Default values for GLMAIN
#
ics = 1  , mxx = 50 , ilg = 1 , ipo = 1 , iug = 1), envir=.GlobalEnv, inherits=TRUE)
}
