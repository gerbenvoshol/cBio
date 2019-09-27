## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libcbio.la

libcbio_la_CFLAGS = $(WARN_CFLAGS) $(DEVWARN_CFLAGS) $(JAVA_CFLAGS) -DPREFIX=\"$(prefix)\" -DBUILD_DIR=\"$(srcdir)\" -DDRV_DIR=\"$(srcdir)\" -DEMBOSS_TOP=\"$(EMBOSS_TOP)\"

libcbio_la_CPPFLAGS = $(JAVA_CPPFLAGS)
libcbio_la_CPPFLAGS += $(MYSQL_CPPFLAGS)
libcbio_la_CPPFLAGS += $(POSTGRESQL_CPPFLAGS)
libcbio_la_CPPFLAGS += -DAJAX_FIXED_ROOT=$(AJAX_FIXED_ROOT)
libcbio_la_CPPFLAGS += -DAJAX_SYSTEM="$(AJAX_SYSTEM)"
libcbio_la_CPPFLAGS += -DPREFIX=\"$(prefix)\"

libcbio_la_LDFLAGS  = -version-info 6:6:0
libcbio_la_LDFLAGS += $(JAVA_LDFLAGS)
libcbio_la_LDFLAGS += $(MYSQL_LDFLAGS)
libcbio_la_LDFLAGS += $(POSTGRESQL_LDFLAGS)
if ISCYGWIN
libcbio_la_LDFLAGS += -no-undefined
libcbio_la_LDFLAGS += -Wl,--add-stdcall-alias
libcbio_la_LDFLAGS += -Wl,--export-all-symbols
libcbio_la_LDFLAGS += -Wl,--enable-auto-import
endif

libcbio_la_SOURCES = \
adler32.c \
ajacd.c \
ajalign.c \
ajarr.c \
ajassem.c \
ajassemdb.c \
ajassemread.c \
ajassemwrite.c \
ajassert.c \
ajbamindex.c \
ajbase.c \
ajcall.c \
ajcigar.c \
ajcod.c \
ajdan.c \
ajdas.c \
ajdmx.c \
ajdomain.c \
ajdom.c \
ajexcept.c \
ajfeat.c \
ajfeatdb.c \
ajfeatread.c \
ajfeatwrite.c \
ajfile.c \
ajfiledata.c \
ajfileio.c \
ajfmt.c \
ajftp.c \
ajgraph.c \
ajhist.c \
ajhttp.c \
ajindex.c \
ajjava.c \
ajlist.c \
ajmart.c \
ajmath.c \
ajmatrices.c \
ajmem.c \
ajmess.c \
ajnam.c \
ajnexus.c \
ajobo.c \
ajobodb.c \
ajoboread.c \
ajobowrite.c \
ajpat.c \
ajpdb.c \
ajpdbio.c \
ajphylo.c \
ajquery.c \
ajrange.c \
ajrefseq.c \
ajrefseqdb.c \
ajrefseqread.c \
ajrefseqwrite.c \
ajreg.c \
ajreport.c \
ajresource.c \
ajresourcedb.c \
ajresourceread.c \
ajresourcewrite.c \
ajseqabi.c \
ajseqbam.c \
ajseq.c \
ajseqdb.c \
ajseqread.c \
ajseqtype.c \
ajseqwrite.c \
ajsoap.c \
ajsort.c \
ajsql.c \
ajstr.c \
ajsys.c \
ajtable.c \
ajtagval.c \
ajtax.c \
ajtaxdb.c \
ajtaxread.c \
ajtaxwrite.c \
ajtext.c \
ajtextdb.c \
ajtextread.c \
ajtextwrite.c \
ajtime.c \
ajtranslate.c \
ajtree.c \
ajurl.c \
ajurldb.c \
ajurlread.c \
ajurlwrite.c \
ajutil.c \
ajvarbcf.c \
ajvar.c \
ajvardb.c \
ajvarread.c \
ajvarwrite.c \
ajvector.c \
ajxml.c \
ajxmldb.c \
ajxmlread.c \
ajxmlwrite.c \
compress.c \
crc32.c \
data.c \
deflate.c \
embaln.c \
embcom.c \
embcons.c \
embdata.c \
embdbi.c \
embdmx.c \
embdomain.c \
embest.c \
embexit.c \
embgroup.c \
embiep.c \
embindex.c \
embinit.c \
embmat.c \
embmisc.c \
embmol.c \
embnmer.c \
embpat.c \
embpatlist.c \
embpdb.c \
embprop.c \
embread.c \
embshow.c \
embsig.c \
embword.c \
ensalign.c \
ensanalysis.c \
ensassembly.c \
ensassemblyexception.c \
ensassemblymapper.c \
ensattribute.c \
ensbaseadaptor.c \
enscache.c \
enscoordsystem.c \
ensdatabaseadaptor.c \
ensdatabaseconnection.c \
ensdatabaseentry.c \
ensdensity.c \
ensditag.c \
ensexon.c \
ensexternaldatabase.c \
ensfeature.c \
ensgene.c \
ensgvallele.c \
ensgvattribute.c \
ensgvbaseadaptor.c \
ensgvdatabaseadaptor.c \
ensgvgenotype.c \
ensgvindividual.c \
ensgvpopulation.c \
ensgvsample.c \
ensgvsource.c \
ensgvsynonym.c \
ensgvvariation.c \
ensintron.c \
enskaryotype.c \
ensmapper.c \
ensmarker.c \
ensmetacoordinate.c \
ensmetainformation.c \
ensmiscellaneous.c \
ensontology.c \
ensoperon.c \
ensprediction.c \
ensprojectionsegment.c \
ensprotein.c \
ensqcalignment.c \
ensqc.c \
ensqcdasfeature.c \
ensqcdatabase.c \
ensqcsequence.c \
ensqcsubmission.c \
ensqcvariation.c \
ensregistry.c \
ensrepeat.c \
ensseqregion.c \
enssequence.c \
enssequenceedit.c \
enssimple.c \
ensslice.c \
enssoftware.c \
ensstorable.c \
enstable.c \
enstranscript.c \
enstranslation.c \
ensutility.c \
ensvariation.c \
gd.c \
gzclose.c \
gzlib.c \
gzread.c \
gzwrite.c \
hpgl.c \
impress.c \
infback.c \
inffast.c \
inflate.c \
inftrees.c \
ljii.c \
ljiip.c \
null.c \
pbm.c \
pcre_chartables.c \
pcre_compile.c \
pcre_dfa_exec.c \
pcre_exec.c \
pcre_fullinfo.c \
pcre_get.c \
pcre_globals.c \
pcre_info.c \
pcre_newline.c \
pcreposix.c \
pcre_study.c \
pcre_tables.c \
pcre_try_flipped.c \
pcre_version.c \
pcre_xclass.c \
pdf.c \
pdfutils.c \
plargs.c \
plbox.c \
plbuf.c \
plcont.c \
plcore.c \
plctrl.c \
plcvt.c \
pldtik.c \
plfill.c \
plhist.c \
plimage.c \
plline.c \
plmap.c \
plmeta.c \
plot3d.c \
plpage.c \
plsdef.c \
plshade.c \
plstdio.c \
plstripc.c \
plsym.c \
pltick.c \
plvect.c \
plvpor.c \
plwind.c \
ps.c \
svg.c \
tek.c \
trees.c \
uncompr.c \
xfig.c \
xmlparse.c \
xmlrole.c \
xmltok.c \
xmltok_impl.c \
xmltok_ns.c \
xwin.c \
zutil.c\
libbwa.c\
libbwa_index.c\
libbwa_aln.c\
libbwa_samse.c\
libbwa_sampe.c\
libbwa_sw.c\
libbwa_mem.c\
libbwa_fastmap.c\
libbwa_fa2pac.c\
libbwa_pac2bwt.c\
libbwa_bwtgen.c\
utils.c\
kthread.c\
kstring.c\
ksw.c\
bwt.c\
bntseq.c\
bwa.c\
bwamem.c\
bwamem_pair.c\
bwamem_extra.c\
malloc_wrap.c\
QSufSort.c\
bwt_gen.c\
bwase.c\
bwaseqio.c\
bwtgap.c\
bwtaln.c\
bamlite.c\
is.c\
bwtindex.c\
bwape.c\
kopen.c\
pemerge.c\
bwtsw2_core.c\
bwtsw2_main.c\
bwtsw2_aux.c\
bwt_lite.c\
bwtsw2_chain.c\
fastmap.c\
bwtsw2_pair.c\
easellib.c\
hmmerlib.c\
infernallib.c


include_HEADERS = \
ajacd.h\
ajalign.h\
ajarch.h\
ajarr.h\
ajassemdata.h\
ajassemdb.h\
ajassem.h\
ajassemread.h\
ajassemwrite.h\
ajassert.h\
ajax.h\
ajbamindex.h\
ajbase.h\
ajcall.h\
ajcigar.h\
ajcod.h\
ajdan.h\
ajdas.h\
ajdefine.h\
ajdmx.h\
ajdomain.h\
ajdom.h\
ajexcept.h\
ajfeatdata.h\
ajfeatdb.h\
ajfeat.h\
ajfeatread.h\
ajfeatwrite.h\
ajfiledata.h\
ajfile.h\
ajfileio.h\
ajfmt.h\
ajftp.h\
ajgraph.h\
ajgraphstruct.h\
ajhist.h\
ajhttp.h\
ajindex.h\
ajjava.h\
ajlib.h\
ajlist.h\
ajmart.h\
ajmath.h\
ajmatrices.h\
ajmem.h\
ajmess.h\
ajnam.h\
ajnexus.h\
ajobodata.h\
ajobodb.h\
ajobo.h\
ajoboread.h\
ajobowrite.h\
ajpat.h\
ajpdb.h\
ajpdbio.h\
ajphylo.h\
ajquerydata.h\
ajquery.h\
ajrange.h\
ajrefseqdata.h\
ajrefseqdb.h\
ajrefseq.h\
ajrefseqread.h\
ajrefseqwrite.h\
ajreg.h\
ajreport.h\
ajresourcedata.h\
ajresourcedb.h\
ajresource.h\
ajresourceread.h\
ajresourcewrite.h\
ajseqabi.h\
ajseqbam.h\
ajseqdata.h\
ajseqdb.h\
ajseq.h\
ajseqread.h\
ajseqtype.h\
ajseqwrite.h\
ajsoap.h\
ajsort.h\
ajsql.h\
ajstr.h\
ajsys.h\
ajtable.h\
ajtagval.h\
ajtaxdata.h\
ajtaxdb.h\
ajtax.h\
ajtaxread.h\
ajtaxwrite.h\
ajtextdata.h\
ajtextdb.h\
ajtext.h\
ajtextread.h\
ajtextwrite.h\
ajtime.h\
ajtranslate.h\
ajtree.h\
ajurldata.h\
ajurldb.h\
ajurl.h\
ajurlread.h\
ajurlwrite.h\
ajutil.h\
ajvarbcf.h\
ajvardata.h\
ajvardb.h\
ajvar.h\
ajvarread.h\
ajvarwrite.h\
ajvector.h\
ajxmldata.h\
ajxmldb.h\
ajxml.h\
ajxmlread.h\
ajxmlwrite.h\
ascii.h\
asciitab.h\
config.h\
crc32.h\
deflate.h\
disptab.h\
drivers.h\
embaln.h\
embcom.h\
embcons.h\
embdata.h\
embdbi.h\
embdmx.h\
embdomain.h\
embest.h\
embexit.h\
embgroup.h\
embiep.h\
embindex.h\
embinit.h\
embmat.h\
embmisc.h\
embmol.h\
embnmer.h\
cbio.h\
embpat.h\
embpatlist.h\
embpdb.h\
embprop.h\
embread.h\
embshow.h\
embsig.h\
embword.h\
ensalign.h\
ensanalysis.h\
ensassemblyexception.h\
ensassembly.h\
ensassemblymapper.h\
ensattribute.h\
ensbaseadaptor.h\
enscache.h\
enscoordsystem.h\
ensdatabaseadaptor.h\
ensdatabaseconnection.h\
ensdatabaseentry.h\
ensdata.h\
ensdensity.h\
ensditag.h\
ensembl.h\
ensexon.h\
ensexternaldatabase.h\
ensfeature.h\
ensgene.h\
ensgvallele.h\
ensgvattribute.h\
ensgvbaseadaptor.h\
ensgvdatabaseadaptor.h\
ensgvdata.h\
ensgvgenotype.h\
ensgvindividual.h\
ensgvpopulation.h\
ensgvsample.h\
ensgvsource.h\
ensgvsynonym.h\
ensgvvariation.h\
ensintron.h\
enskaryotype.h\
ensmapper.h\
ensmarker.h\
ensmetacoordinate.h\
ensmetainformation.h\
ensmiscellaneous.h\
ensontology.h\
ensoperon.h\
ensprediction.h\
ensprojectionsegment.h\
ensprotein.h\
ensqcalignment.h\
ensqcdasfeature.h\
ensqcdatabase.h\
ensqc.h\
ensqcsequence.h\
ensqcsubmission.h\
ensqcvariation.h\
ensregistry.h\
ensrepeat.h\
ensseqregion.h\
enssequenceedit.h\
enssequence.h\
enssimple.h\
ensslice.h\
enssoftware.h\
ensstorable.h\
enstable.h\
enstranscript.h\
enstranslation.h\
ensutility.h\
ensvariation.h\
expat_config.h\
expatemboss.h\
expat_external.h\
expat.h\
gzguts.h\
iasciitab.h\
inffast.h\
inffixed.h\
inflate.h\
inftrees.h\
internal.h\
latin1tab.h\
metadefs.h\
nametab.h\
pcre_config.h\
pcreemboss.h\
pcre.h\
pcre_internal.h\
pcreposix.h\
pdf.h\
plConfig.h\
plcore.h\
pldebug.h\
plDevs.h\
pldll.h\
plevent.h\
plfci-type1.h\
plhershey-unicode.h\
plplot.h\
plplotP.h\
plstrm.h\
plunicode-type1.h\
plxwd.h\
ps.h\
trees.h\
ucp.h\
utf8tab.h\
xmlrole.h\
xmltok.h\
xmltok_impl.h\
zconf.h\
zlibemboss.h\
zlib.h\
zutil.h\
libbwa.h