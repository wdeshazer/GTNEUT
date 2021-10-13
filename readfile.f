      program readfile
        implicit none

        include 'neutGlob.inc'
        include 'comiou.inc'
        include 'consts.inc'

        namelist /inp/ i_inp, nCells, nWallSegm, nPlasmReg, iType,
     .   nSides, lside, angle, adjCell, scalFact, idbug,
     .   elecTemp, ionTemp, elecDens, ionDens, S_ext, Rwall, g_ion,
     .   i_e0, eneut, eneut_v, zion, aion, aneut, g_ex, prntOrdr,
     .   iatdat, lchex, leh0, v0fact, icosn, iquad, nph,
     .   iescp, irefl, ifjsv, twall, fwabsorb, zwall, awall, ifrstcol,
     .   isparsitr,idp, inon, nd0, neitr , Shotnumber, Timeslice,nxleg1,
     .   nxcore1,nxcore2,nycore1,nysol1,nxxpt,nxmod,nxleg2
 
 
        namelist /inp1/ NX, NY, Lx, Ly, te_fixed, ti_fixed, ne_fixed,
     .   ni_fixed, S_0, r_lft, r_rgt, r_top, r_btm, g_lft, g_rgt, g_top,
     .   g_btm, flx_lft, flx_rgt, flx_top, flx_btm, igradneh, igradnev,
     .   igradnih, igradniv, igradteh, igradtev, igradtih, igradtiv,
     .   ne_lft, ne_rgt, ne_top, ne_btm, ni_lft, ni_rgt, ni_top, ni_btm,
     .   te_lft, te_rgt, te_top, te_btm, ti_lft, ti_rgt, ti_top, ti_btm,
     .   iexp, ialphaAll, ialphaDen, alpha

        open (nin, file = 'toneut', status = 'old')
   
        read (nin,inp)
   
        if (i_inp.EQ.1) read (nin, inp1)
   
      stop
      end