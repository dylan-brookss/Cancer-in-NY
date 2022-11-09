#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2020
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp.
# ************************************************************************/"""Partial Least Squares Regression Module"""



import random, os, tempfile, textwrap

"""STATS CANCORR extension command"""

__author__ =  'IBM SPSS, JKP'
__version__=  '1.0.1'

# history
# 04-jun-2013 Original version

helptext = """STATS CANCORR SET1=variable list SET2=variable list
/OPTIONS
    SCORINGSYNTAX = filespec
    ROOTNAME = variable name
    COMPUTECVARS =YES or NO*
    CENTEREDDS = dataset name
    NDIMS = integer
    
/PRINT PAIRWISECORR=YES* or NO
    LOADINGS = YES* or NO
    VARPROP = YES* or NO
    COEFFICIENTS = YES* or NO
[/HELP]
* indicates default value.

This command computes canonical correlations between two sets
of variables and optionally create a dataset of scores.  Listwise
deletion is used for missing values.

SET1 and SET2 specify the two sets of variables whose
canonical correlations should be computed.

If SCORINGSYNTAX is specified, a syntax file creating code to
calcuate scores is written.  It can be used with another dataset
if the names match.  The syntax begins with mean centering the
input data. ROOTNAME is required if this keyword is used.

ROOTNAME specifies a prefix for the score variables to be created.
The full names are ROOTNAME_set#_corr# where set# is set1 or set2 and
corr# goes from 1 to the number of canonical correlations.  Keep
the root name short enough that the complete name is no longer
than 64 bytes.

COMPUTECVARS specifies whether or not to compute the scores.
If YES, the input variables are mean centered and written to a
new dataset which also contains the scores.  CENTEREDDS and
ROOTNAME are required if this is used.

CENTEREDDS is the name for the centered dataset.  It is
ignored if COMPUTECVARS is NO.

By default, all dimensions of the canonical correlations are
used for scoring.  Specify NDIMS=n to use a smaller number.


The PRINT subcommand controls what optional output is
displayed.

Example:
STATS CANCOR SET1=y1 y2 SET2=x1 x2 x3
/OPTIONS SCORINGSYNTAX = "C:/syntax/scoring.sps"
ROOTNAME = score COMPUTECVARS = YES CENTEREDDS = centered.

/HELP displays this help and does nothing else.
"""

import spss, spssaux, spssdata
from extension import Template, Syntax, processcmd

def getTempDir(tempdir=None):
    """Return a temporary directory"""
    
    
    trythese = ["SPSSTEMP", "TEMP", "TMP", "TMPDIR"]
    if tempdir is None:
        for d in trythese:
            try:
                tname = os.environ[d]
                break
            except:
                continue
        else:
            tname = tempfile.gettempdir()
    else:
        tname = tempdir
    tname = tname.replace("\\", "/")
    if tname[-1] == "/":   # remove any trailing separator
        tname = tname[:-1]
    return tname

def mkrandomname(prefix="D", sav=True):
    res = prefix + str(random.random())
    if sav:
        res = res + ".sav"
    return res

def cancorr(set1, set2, pairwisecorr=True, loadings=True, 
    varprop=True, coefficients=True,
    computecvars = False, scoringsyntax=None, rootname=None, 
    centeredds=None, ndims=None):
    """Compute and display canonical correlations"""
    
    # debugging
    # makes debug apply only to the current thread
    #try:
        #import wingdbstub
        #if wingdbstub.debugger != None:
            #import time
            #wingdbstub.debugger.StopDebug()
            #time.sleep(1)
            #wingdbstub.debugger.StartDebug()
        #import thread
        #wingdbstub.debugger.SetDebugThreads({thread.get_ident(): 1}, default_policy=0)
        ## for V19 use
        ###    ###SpssClient._heartBeat(False)
    #except:
        #pass
    
    # check variable types
    types = spssaux.VariableDict(set1+set2, variableType="string")
    if len(types) > 0:   # VariableDict objects implement len as the number of variables in the dictionary
        raise ValueError(_("String variables are not allowed in this procedure"))
    s = set(set1).intersection(set2)
    if s:
        raise ValueError(_("""The two sets of variables have variables in common: %s""") % " ".join(s))
    if (scoringsyntax or computecvars) and not rootname:
        raise ValueError(_("A root name for output variables is required if computing scores or saving syntax"))
    if ((scoringsyntax or computecvars) and centeredds is None):
        raise ValueError(_("""A name for the centered dataset is required if scoring or creating scoring syntax"""))
    
    activeds = spss.ActiveDataset()
    tempdir = getTempDir()
    if activeds == "*":
        activeds = mkrandomname(sav=False)
        spss.Submit("DATASET NAME %s" % activeds)
    fh = spssaux.FileHandles()
    if scoringsyntax:
        fh = fh.resolve(scoringsyntax)
    resdict = docalculations(set1, set2, tempdir, pairwisecorr, activeds)
    displayresults(set1, set2, resdict, centeredds, scoringsyntax,
        loadings, varprop, coefficients, ndims)

    # Calculate scores and write scoring syntax if indicated
    scoresandsyntax(scoringsyntax, computecvars, rootname, centeredds, resdict, 
        set1, set2, ndims)


def scoresandsyntax(scoringsyntax, computecvars, rootname, centeredds, resdict, 
        set1, set2, ndims):
    """Calculate scores and write scoring syntax file if specified"""
    
    #a1 is unstandardized c coefs, set1
    #b1 is unstandardized c coefs, set2
    
    if not (computecvars or scoringsyntax):
        return
    lhsvars = set()
    syntaxlist = []
    syntax, newnames = buildsyntax(rootname, 1, set1, resdict['a1'], ndims)
    lhsvars.update(newnames)
    syntaxlist.append(syntax)
    syntax, newnames = buildsyntax(rootname, 2, set2, resdict['b1'], ndims)
    lhsvars.update(newnames)
    syntaxlist.append(syntax)
    if lhsvars.intersection(set(set1+set2)):
        raise ValueError(_("""A scoring variable name appears in the list of input variables"""))
    centeringsyntax = r"""DATASET DECLARE %(centeredds)s.
* Center the variables.
matrix.
get M /variables %(variables)s/missing=omit.
compute means = csum(M) / NROW(M).
compute numcol = NCOL(M).
loop i = 1 to numcol.
compute M(:,i) = M(:,i) - means(i).
end loop.
save M /outfile=%(centeredds)s /variables=%(variables)s.
end matrix.
dataset activate %(centeredds)s.
""" % {"variables" : " ".join(set1 + set2), "centeredds" : centeredds}

    if scoringsyntax:
        f = open(scoringsyntax, "w")
        f.write(_("""* Calculate canonical scores based on canonical correlations.\n"""))
        f.write(centeringsyntax)
        f.write("\n".join(syntaxlist))  # does not add newline characters
        f.close()
        
    if computecvars:
        spss.EndProcedure()
        # leave calculations exposed in case of errors
        spss.Submit(centeringsyntax)
        spss.Submit(syntaxlist)

def docalculations(set1, set2, tempdir, pairwisecorr, activeds):
    """compute canonical correlations and return results dictionary
    
    tempdir is where to write the temporary sav files
    printcorr is whether or not to print correlations
    activeds is the name of the active dataset"""
    
    # In order to create pivot table output from MATRIX results
    # we write a temporary sav file from MATRIX for each ouput and erase it
    # later
    
    set1str = " ".join(set1)
    set2str = " ".join(set2)    

    test = mkrandomname(tempdir + "/")
    a = mkrandomname(tempdir + "/")  #set1 cc
    a1 = mkrandomname(tempdir + "/")  #set1 cc standardized
    b = mkrandomname(tempdir + "/")  #set2 cc
    b1 = mkrandomname(tempdir + "/") #set2 cc standardized
    f1 = mkrandomname(tempdir + "/") #redundancy index
    tem = mkrandomname(tempdir + "/")
    tem2 = mkrandomname(tempdir + "/") # cross loadings set1
    tem3 = mkrandomname(tempdir + "/") # cross loadings set2
    tem4 = mkrandomname(tempdir + "/") # 
    redun = mkrandomname(tempdir + "/") # f1, cs3, f2, cs4
    canonical_correlation_size = mkrandomname(tempdir + "/")
    canonical_correlation_AB = mkrandomname(tempdir + "/")
    corrmatfile = mkrandomname(tempdir + "/")
    omsrequest = mkrandomname("O")
    if not pairwisecorr:
        spss.Submit(r"""OMS SELECT ALL EXCEPT=WARNINGS
        /DESTINATION VIEWER=NO /TAG = %(omsrequest)s.""" % locals())
    spss.Submit(r"""CORRELATIONS %(set1str)s %(set2str)s /MISSING=LISTWISE/PRINT=SIG
    /MATRIX OUT('%(corrmatfile)s').""" % locals())

    # If output not already suppressed, suppress most of it now.
    if pairwisecorr:
        spss.Submit(r"""OMS SELECT ALL EXCEPT=[WARNINGS]
        /DESTINATION VIEWER=NO /TAG = %(omsrequest)s.""" % locals())
    try:
        # We can't continue a MATRIX command across Submit calls 
        spss.Submit("""PRESERVE.
SET MXLOOPS = 500 MITERATE 500.
matrix.
get r /file = "%(corrmatfile)s" /variables = %(set1str)s.
compute p1 = ncol(r).
get r /file = "%(corrmatfile)s" /names = varname/variables = %(set1str)s %(set2str)s.
compute p2 = ncol(r)-p1.
compute nx1 = varname(1:p1).
compute nv = p1+p2.
compute nx2 = varname((p1+1):nv).
compute rr = r(4:(nv+3),1:nv).
compute ncases = r(3,3).
compute ns = r(3,1).
compute r11 = rr(1:p1,1:p1).
compute r22 = rr((p1+1):nv,(p1+1):nv).
compute r12 = rr(1:p1,(p1+1):nv).
compute d1 = r(2,1:p1).
compute d2 = r(2,(p1+1):nv).
do if (ncol(d1) > 1).
compute d1 = mdiag(d1).
end if.
do if (ncol(d2) > 1).
compute d2 = mdiag(d2).
end if.
compute s1 = d1 * r11 * d1.
compute s12 = d1 * r12 * d2.
compute s2 = d2 * r22 * d2.
compute d1 = inv(d1).
compute d2 = inv(d2).
compute r1 = chol(r11).
compute r2 = chol(r22).
compute r1_inv = inv(r1).
compute r2_inv = inv(r2).
do if (p1 le p2).
compute omega = t(r1_inv) * r12 * r2_inv.
else.
compute omega = t(r2_inv) * t(r12) * r1_inv.
end if.
call svd(omega,u,lambda,v).
compute dlam = diag(lambda).
*  eigenvalues.
compute eign = (1 &/ (1-dlam &** 2)) - 1.
compute wlam = 1 &/ (1+eign).
compute n = nrow(wlam).
compute wilk = wlam.
compute df1 = wlam.
compute df2 = wlam.
compute sig = wlam.
compute bart2 = wlam.
compute temp = 1.
compute p = p1.
compute q = p2.
compute mmm = ncases - 3/2 - (p + q) / 2.
loop i = 1 to n.
   compute temp = temp * wlam(n-i+1).
   compute wilk(n-i+1) = temp.
end loop.
* calculate asymptotic F statistics.
loop  i = 1 to n.
   compute num = p**2 * q**2 - 4.
   do if (num ne 0).
   compute s = sqrt(num / (p**2 + q**2 - 5)).
   compute si = 1/s.
   compute df1(i) = p * q.
   compute df2(i) = mmm * s - p * q/2 + 1.
   compute r = (1 - wilk(i)**si) / wilk(i)**si.
   compute bart2(i) = r * df2(i) / df1(i).
   compute sig(i) = 1 - fcdf(bart2(i), df1(i), df2(i)).
   compute p = p - 1.
   compute q = q - 1.
   else.
   * Can't always compute the F.  This value will be
   * converted to proper sysmis later, since MATRIX
   * does not support sysmis directly.
   compute df1(i) = -99.
   compute df2(i) = -99.
   compute sig(i) = -99.
   compute bart2(i) = -99.
   end if.
+  *compute tem = tem * wlam(n-#l+1).
+  *compute df(n-#l+1) = (p1-n+#l) * (p2-n+#l).
+  *compute dof = df(n-#l+1).
+  *compute bart2(n-#l+1) = -(ns-0.5 * (p1+p2+3)) * ln(tem).
+  *compute chi = bart2(n-#l+1).
+  *compute sig(n-#l+1) = 1-chicdf(chi,dof).
+  *compute wilk(n-i+1) = tem.
end loop.
compute test = {dlam, eign, wilk, bart2, df1, df2, sig}.
save test / outfile = "%(test)s".
*  standardized CCs.
do if (p1 le p2).
compute a = r1_inv * u.
else.
compute a = r1_inv * v.
end if.
do if (p2 lt p1).
compute a = a(:,1:p2).
end if.
*  standardize c coef set1 as a.
save a /outfile = "%(a)s".
*  unstandardized.
compute a1 = d1 * a.
*  unstandardized c coef set1 as a1.
save a1 /outfile = "%(a1)s".
* set2.
do if (p1 le p2).
compute b = r2_inv * v.
else.
compute b = r2_inv * u.
end if.
do if (p1 le p2).
compute b = b(:,1:p1).
end if.
*  unstandardized c coef set2 as b.
save b/outfile = "%(b)s".
compute b1 = d2 * b.
save b1 /outfile = "%(b1)s".
*  loadings -set1.
compute tem = d1 * s1 * a1.
save tem /outfile = "%(tem)s".
*  redundancy index.
compute f1 = cssq(tem)/p1.
compute f1 = t(f1).
save f1 /outfile = "%(f1)s".
*  cross loadings set1 as tem2.
compute tem = d1 * s12 * b1.
save tem /outfile = "%(tem2)s".
*  redundancy index.
compute cs3 = cssq(tem)/p1.
compute cs3 = t(cs3).
compute tem = d2 * s2 * b1.
* canonical loadings set2 as tem4.
save tem /outfile = "%(tem4)s".
compute f2 = cssq(tem)/p2.
compute f2 = t(f2).
*  set2 cross loadings as tem3.
compute tem = d2 * t(s12) * a1.
save tem /outfile = "%(tem3)s".
* redundancy index.
compute cs4 = cssq(tem)/p2.
compute cs4 = t(cs4).
save {f1,cs3,f2,cs4} /outfile = "%(redun)s".
* scoring files.
SAVE {P1,P2} / OUTFILE =  "%(canonical_correlation_size)s".
SAVE {T(A1),T(B1)} / OUTFILE =  "%(canonical_correlation_AB)s".
END MATRIX.""" % locals())

        # Read all results and erase the files
        resdict = dict()
        w = WorkingDatasets()
        reslist = [test, a, a1, b, b1, f1, tem, tem2, tem3, tem4, redun]  #these values are file names
        rn = ['test','a','a1','b','b1','f1','tem','tem2','tem3','tem4','redun']
        for name, f in zip(rn, reslist):
            resdict[name] = w.getsav(f)
    except:
        # Note: this error message will not be seen if SPSS_EXTENSIONS_RAISE is True because of OMS suppression.
        raise ValueError(_("""The requested canonical correlations could not be computed due to singularity or missing data."""))
    finally:
        spss.Submit("RESTORE.")
        spss.Submit("""OMSEND TAG = %s""" % omsrequest)
        spss.Submit("""DATASET ACTIVATE %s.""" % activeds)
    return resdict

def displayresults(set1, set2, resdict, centeredds, scoringsyntax,
        loadings, varprop, coefficients, ndims):
    StartProcedure("Canonical Correlations", "STATSCANCORR")
    ncor = min(len(set1), len(set2))
    if not ndims is None:
        ncorused = min(ndims, ncor)
    else:
        ncorused = ncor
    rowlabels = [_("""Set 1 Variables"""), _("""Set 2 Variables"""), _("""Centered Dataset"""),
        _("""Scoring Syntax"""), _("""Correlations Used for Scoring""")]
    cells = [" ".join(set1), " ".join(set2),
        (centeredds or _("""None""")),
        (scoringsyntax or _("""None""")),
        str(ncorused)]
    t = spss.BasePivotTable(_("Canonical Correlations Settings"), "CANCORRSET")
    t.SimplePivotTable(rowlabels = rowlabels, collabels=[_("Values")],
        cells=cells)
    # canonical correlations
    t = spss.BasePivotTable(_("Canonical Correlations"), "CANCORRCORRS",
            caption=_("""H0 for Wilks test is that the correlations in the current and following rows are zero"""))
    # convert pseudo sysmis value from Matrix to proper sysmis
    # The cases are tuples, so we have to convert to lists to modify.
    resdict['test'] = list(resdict['test'])
    for i in range(ncor):
        resdict['test'][i] = list(resdict['test'][i])
        for j in range(3, 7):
            if resdict['test'][i][j] == -99:
                resdict['test'][i][j] = None

    t.SimplePivotTable(
        rowlabels=[str(i+1) for i in range(ncor)],
        collabels = [_("Correlation"), _("Eigenvalue"), 
            _("Wilks Statistic"), _("F"), _("Num D.F"), _("Denom D.F."), 
            _("Sig.")],
        cells = resdict['test'])

    if coefficients:
        t = spss.BasePivotTable(_("Set 1 Standardized Canonical Correlation Coefficients"), "CANCORRSET1SCORRS")
        t.SimplePivotTable(rowdim = _("Variable"), 
            rowlabels = set1,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['a'])
        # set 2
        t = spss.BasePivotTable(_("Set 2 Standardized Canonical Correlation Coefficients"), "CANCORRSET2SCORRS")
        t.SimplePivotTable(rowdim = _("Variable"), 
            rowlabels = set2,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['b'])
        
        t = spss.BasePivotTable(_("Set 1 Unstandardized Canonical Correlation Coefficients"), "CANCORRSET1UCORRS")
        t.SimplePivotTable(rowdim = _("Variable"), 
            rowlabels = set1,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['a1'])
        t = spss.BasePivotTable(_("Set 2 Unstandardized Canonical Correlation Coefficients"), "CANCORRSET2UCORRS")
        t.SimplePivotTable(rowdim = _("Variable"), 
            rowlabels = set2,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['b1'])
    if loadings:
        t = spss.BasePivotTable(_("Set 1 Canonical Loadings"), "CANCORRSET1LOADINGS")
        t.SimplePivotTable(
            rowdim = _("Variable"),
            rowlabels = set1,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['tem'])
        t = spss.BasePivotTable(_("Set 2 Canonical Loadings"), "CANCORRSET2LOADINGS")
        t.SimplePivotTable(
            rowdim = _("Variable"),
            rowlabels = set2,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['tem4'])
    
        t = spss.BasePivotTable(_("Set 1 Cross Loadings"), "CANCORRSET1CROSSLOADINGS")
        t.SimplePivotTable(
            rowdim = _("Variable"),
            rowlabels = set1,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['tem2'])
        t = spss.BasePivotTable(_("Set 2 Cross Loadings"), "CANCORRSET2CROSSLOADINGS")
        t.SimplePivotTable(
            rowdim = _("Variable"),
            rowlabels = set2,
            collabels = [str(i+1) for i in range(ncor)],
            cells = resdict['tem3'])
    if varprop:
        t = spss.BasePivotTable(_("Proportion of Variance Explained"), "CANCORRREDUN")
        t.SimplePivotTable(
            rowdim=_("Canonical Variable"),
            rowlabels = [str(i+1) for i in range(ncor)],
            collabels=[_("Set 1 by Self"), _("Set 1 by Set 2"), 
                _("Set 2 by Self"), _("Set 2 by Set 1")],
            cells=resdict["redun"])
    
    
class WorkingDatasets(object):
    def __init__(self):
        self.wdsname = mkrandomname("D", sav=False)
        
    def getsav(self, filespec, delete=True):
        """Open sav file and return all contents
        
        filespec is the file path
        filespec is deleted after the contents are read unless delete==False"""
     
        item = self.wdsname
        spss.Submit(r"""get file="%(filespec)s".
DATASET NAME %(item)s.
DATASET ACTIVATE %(item)s.""" % locals())
        contents = spssdata.Spssdata(names=False).fetchall()
        spss.Submit("""DATASET CLOSE %(item)s.
        NEW FILE.""" % locals())
        if delete:
            os.remove(filespec)
        return contents

def buildsyntax(root, setn, setvars, data, ndims):
    """return list of syntax specs
    
    root is the root of the name for the left hand variable
    setn is the set number
    setvars is the list of variables in the set
    data is the table of unstandardized canconical coefficients
        - one column per canonical correlation
    ndims can trim the number of correlations used."""
    
    syntax = []
    nvars = len(setvars)
    ncor = len(data[0])
    if not ndims is None:
        ncor = min(ncor, ndims)
    newnames = set()
    for i in range(ncor):
        cname = root + "_set" + str(setn) + "_" + str(i+1)
        newnames.add(cname)
        if len(cname) > 64:
            raise ValueError(_("The specified root name is too long: %s") % root)
        s = ["COMPUTE " + cname + " = "]
        for j in range(nvars):
            s.append(str(data[j][i]) + " * " + setvars[j])
        syntax.append(s[0] + " + ".join(s[1:]))
        syntax[i] = "\n".join(textwrap.wrap(syntax[i])) +"."

    return "\n".join(syntax), newnames


def Run(args):
    """Execute the STATS CANCORR extension command"""

    args = args[list(args.keys())[0]]

    oobj = Syntax([
        Template("SET1", subc="",  ktype="existingvarlist", var="set1", islist=True),
        Template("SET2", subc="",  ktype="existingvarlist", var="set2", islist=True),

        Template("SCORINGSYNTAX", subc="OPTIONS", ktype="literal", var="scoringsyntax"),
        Template("SCORES", subc="OPTIONS", ktype="literal", var="scores"),
        Template("COMPUTECVARS", subc="OPTIONS", ktype="bool", var="computecvars"),
        Template("ROOTNAME", subc="OPTIONS", ktype="literal", var="rootname"),
        Template("CENTEREDDS", subc="OPTIONS", ktype="varname", var="centeredds"),
        Template("NDIMS", subc="OPTIONS", ktype="int", var="ndims", vallist=(1,)),
        
        Template("PAIRWISECORR", subc="PRINT", ktype="bool", var="pairwisecorr"),
        Template("LOADINGS", subc="PRINT", ktype="bool", var="loadings"),
        Template("VARPROP", subc="PRINT", ktype="bool", var="varprop"),
        Template("COEFFICIENTS", subc="PRINT", ktype="bool", var="coefficients"),
        
        Template("HELP", subc="", ktype="bool")])
    
    #enable localization
    global _
    try:
        _("---")
    except:
        def _(msg):
            return msg
    # A HELP subcommand overrides all else
    if "HELP" in args:
        #print helptext
        helper()
    else:
        processcmd(oobj, args, cancorr, vardict=spssaux.VariableDict())

def helper():
    """open html help in default browser window
    
    The location is computed from the current module name"""
    
    import webbrowser, os.path
    
    path = os.path.splitext(__file__)[0]
    helpspec = "file://" + path + os.path.sep + \
         "markdown.html"
    
    # webbrowser.open seems not to work well
    browser = webbrowser.get()
    if not browser.open_new(helpspec):
        print(("Help file not found:" + helpspec))
try:    #override
    from extension import helper
except:
    pass        
class NonProcPivotTable(object):
    """Accumulate an object that can be turned into a basic pivot table once a procedure state can be established"""
    
    def __init__(self, omssubtype, outlinetitle="", tabletitle="", caption="", rowdim="", coldim="", columnlabels=[],
                 procname="Messages"):
        """omssubtype is the OMS table subtype.
        caption is the table caption.
        tabletitle is the table title.
        columnlabels is a sequence of column labels.
        If columnlabels is empty, this is treated as a one-column table, and the rowlabels are used as the values with
        the label column hidden
        
        procname is the procedure name.  It must not be translated."""
        
        attributesFromDict(locals())
        self.rowlabels = []
        self.columnvalues = []
        self.rowcount = 0

    def addrow(self, rowlabel=None, cvalues=None):
        """Append a row labelled rowlabel to the table and set value(s) from cvalues.
        
        rowlabel is a label for the stub.
        cvalues is a sequence of values with the same number of values are there are columns in the table."""

        if cvalues is None:
            cvalues = []
        self.rowcount += 1
        if rowlabel is None:
            self.rowlabels.append(str(self.rowcount))
        else:
            self.rowlabels.append(rowlabel)
        self.columnvalues.extend(cvalues)
        
    def generate(self):
        """Produce the table assuming that a procedure state is now in effect if it has any rows."""
        
        privateproc = False
        if self.rowcount > 0:
            try:
                table = spss.BasePivotTable(self.tabletitle, self.omssubtype)
            except:
                StartProcedure(_("Create dummy variables"), self.procname)
                privateproc = True
                table = spss.BasePivotTable(self.tabletitle, self.omssubtype)
            if self.caption:
                table.Caption(self.caption)
            if self.columnlabels != []:
                table.SimplePivotTable(self.rowdim, self.rowlabels, self.coldim, self.columnlabels, self.columnvalues)
            else:
                table.Append(spss.Dimension.Place.row,"rowdim",hideName=True,hideLabels=True)
                table.Append(spss.Dimension.Place.column,"coldim",hideName=True,hideLabels=True)
                colcat = spss.CellText.String("Message")
                for r in self.rowlabels:
                    cellr = spss.CellText.String(r)
                    table[(cellr, colcat)] = cellr
            if privateproc:
                spss.EndProcedure()
                
def attributesFromDict(d):
    """build self attributes from a dictionary d."""
    self = d.pop('self')
    for name, value in d.items():
        setattr(self, name, value)

def StartProcedure(procname, omsid):
    """Start a procedure
    
    procname is the name that will appear in the Viewer outline.  It may be translated
    omsid is the OMS procedure identifier and should not be translated.
    
    Statistics versions prior to 19 support only a single term used for both purposes.
    For those versions, the omsid will be use for the procedure name.
    
    While the spss.StartProcedure function accepts the one argument, this function
    requires both."""
    
    try:
        spss.StartProcedure(procname, omsid)
    except TypeError:  #older version
        spss.StartProcedure(omsid)