"dfcomn" <-
function(ipsi=-9,c=-1.345,h1=-1.7,h2=-3.4,h3=-8.5,xk=-1.548,d=-1.345,
                   beta=-0.5,bet0=-1.0,iucv=-1,a2=0.,b2=-3.0,chk=-9.0,ckw=-2.0,
                   bb=-1.0,bt=-1.0,cw=-1.0,em=-1.345,cr=-2.0,vk=-1.0,np=-2,
                   nu=-1, v7=-1, iwww=-1)
{
        if (!exists(".dFv", , envir=.GlobalEnv)) dfvals()
        f.res <- .Fortran("dfcomn",
        ipsi = to.integer(ipsi),
        c = to.single(c),
        h1 = to.single(h1),
        h2 = to.single(h2),
        h3 = to.single(h3),
        xk = to.single(xk),
        d = to.single(d),
        beta = to.single(beta),
        bet0 = to.single(bet0),
        iucv = to.integer(iucv),
        a2 = to.single(a2),
        b2 = to.single(b2),
        chk = to.single(chk),
        ckw = to.single(ckw),
        bb = to.single(bb),
        bt = to.single(bt),
        cw = to.single(cw),
        em = to.single(em),
        cr = to.single(cr),
        vk = to.single(vk),
        np = to.integer(np),
        enu = to.single(nu),
        v7 = to.single(v7),
        iwww = to.integer(iwww))
        .def <- get(".dFv", envir=.GlobalEnv)
        if (c > 0) .def$ccc <- c 
        if (.def$ccc < 0) .def$ccc <- 1.345
        if (d > 0) .def$ddd <- d 
        if (.def$ddd < 0) .def$ddd <- .def$ccc
#       .dFv <<- .def
        assign(".dFv", .def, envir=.GlobalEnv)
        list(ipsi=f.res$ipsi,iucv=f.res$iucv,iwww=f.res$iwww)
}

