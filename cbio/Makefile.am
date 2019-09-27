## Process this file with automake to produce Makefile.in

SUBDIRS = acd data index

AM_CFLAGS  = $(WARN_CFLAGS) $(DEVWARN_CFLAGS)
if ISSHARED
if ISAIXIA64
AM_CFLAGS += -Wl,-bdynamic
AM_CFLAGS += -Wl,-brtl
AM_CFLAGS += -L../libcbio/.libs
AM_CFLAGS += -lcbio -lm -lgd
endif
endif

AM_CPPFLAGS  = -I$(top_srcdir)/libcbio
AM_CPPFLAGS += $(AXIS2C_CPPFLAGS)
AM_CPPFLAGS += $(JAVA_CPPFLAGS)

LDADD  = ../libcbio/libcbio.la
LDADD += $(XLIB)

#if PURIFY
#LINK = purify $(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(LDFLAGS) -o $@
#endif

# Applications, which are not built.
#
# dbget
# dbgeturl
# dbshowdat
# dbshowquery
# dbtellquery
# docxmlvalid
# docxmlwell
# entrygo
# entrytax
# goseq
# isdbdata
# isdbname
# seqgo
# seqtax
# taxseq
# template ###
# texttotext
# texttoxml
# textvalid
# xmltotext
# xmltoxml
# xmltransform

check_PROGRAMS = \
cB_acdrelations \
cB_acdxsd \
cB_ajtest \
cB_complex \
cB_corbatest \
cB_dastest \
cB_dbfetchtest \
cB_domtesta \
cB_domtestb \
cB_domtestc \
cB_domtestd \
cB_edamclean \
cB_ensembltest \
cB_entrails \
cB_entrailsbook \
cB_entrailshtml \
cB_entrailswiki \
cB_histogramtest \
cB_intconv \
cB_martattributes \
cB_martdatasets \
cB_martfilters \
cB_martquery \
cB_martregistry \
cB_martseqs \
cB_newcoils \
cB_origsplitter \
cB_origunion \
cB_patmattest \
cB_prima \
cB_primers \
cB_seqinfo \
cB_seqretall \
cB_seqretallfeat \
cB_seqretset \
cB_seqretsingle \
cB_seqrettype \
cB_sqltest \
cB_testplot \
cB_treetypedisplay

bin_PROGRAMS = \
cB_aaindexextract \
cB_abiview \
cB_acdc \
cB_acdgalaxy \
cB_acdlog \
cB_acdpretty \
cB_acdtable \
cB_acdtrace \
cB_acdvalid \
cB_aligncopy \
cB_aligncopypair \
cB_antigenic \
cB_assemblyget \
cB_backtranambig \
cB_backtranseq \
cB_banana \
cB_biosed \
cB_btwisted \
cB_cachedas \
cB_cachedbfetch \
cB_cacheebeyesearch \
cB_cacheensembl \
cB_cai \
cB_chaos \
cB_charge \
cB_checktrans \
cB_chips \
cB_cirdna \
cB_codcmp \
cB_codcopy \
cB_coderet \
cB_compseq \
cB_cons \
cB_consambig \
cB_cpgplot \
cB_cpgreport \
cB_cusp \
cB_cutgextract \
cB_cutseq \
cB_dan \
cB_dbiblast \
cB_dbifasta \
cB_dbiflat \
cB_dbigcg \
cB_dbtell \
cB_dbxcompress \
cB_dbxedam \
cB_dbxfasta \
cB_dbxflat \
cB_dbxgcg \
cB_dbxobo \
cB_dbxreport \
cB_dbxresource \
cB_dbxstat \
cB_dbxtax \
cB_dbxuncompress \
cB_degapseq \
cB_density \
cB_descseq \
cB_diffseq \
cB_distmat \
cB_dotmatcher \
cB_dotpath \
cB_dottup \
cB_dreg \
cB_drfinddata \
cB_drfindformat \
cB_drfindid \
cB_drfindresource \
cB_drget \
cB_drtext \
cB_edamdef \
cB_edamhasinput \
cB_edamhasoutput \
cB_edamisformat \
cB_edamisid \
cB_edamname \
cB_edialign \
cB_einverted \
cB_embossdata \
cB_embossupdate \
cB_embossversion \
cB_emma \
cB_emowse \
cB_entret \
cB_epestfind \
cB_eprimer3 \
cB_eprimer32 \
cB_equicktandem \
cB_est2genome \
cB_etandem \
cB_extractalign \
cB_extractfeat \
cB_extractseq \
cB_featcopy \
cB_featmerge \
cB_featreport \
cB_feattext \
cB_findkm \
cB_freak \
cB_fuzznuc \
cB_fuzzpro \
cB_fuzztran \
cB_garnier \
cB_geecee \
cB_getorf \
cB_godef \
cB_goname \
cB_helixturnhelix \
cB_hmoment \
cB_iep \
cB_infoalign \
cB_infoassembly \
cB_infobase \
cB_inforesidue \
cB_infoseq \
cB_isochore \
cB_jaspextract \
cB_jaspscan \
cB_jembossctl \
cB_lindna \
cB_listor \
cB_makenucseq \
cB_makeprotseq \
cB_marscan \
cB_maskambignuc \
cB_maskambigprot \
cB_maskfeat \
cB_maskseq \
cB_matcher \
cB_megamerger \
cB_merger \
cB_msbar \
cB_mwcontam \
cB_mwfilter \
cB_needle \
cB_needleall \
cB_newcpgreport \
cB_newcpgseek \
cB_newseq \
cB_nohtml \
cB_noreturn \
cB_nospace \
cB_notab \
cB_notseq \
cB_nthseq \
cB_nthseqset \
cB_octanol \
cB_oddcomp \
cB_ontocount \
cB_ontoget \
cB_ontogetcommon \
cB_ontogetdown \
cB_ontogetobsolete \
cB_ontogetroot \
cB_ontogetsibs \
cB_ontogetup \
cB_ontoisobsolete \
cB_ontotext \
cB_palindrome \
cB_pasteseq \
cB_patmatdb \
cB_patmatmotifs \
cB_pepcoil \
cB_pepdigest \
cB_pepinfo \
cB_pepnet \
cB_pepstats \
cB_pepwheel \
cB_pepwindow \
cB_pepwindowall \
cB_plotcon \
cB_plotorf \
cB_polydot \
cB_preg \
cB_prettyplot \
cB_prettyseq \
cB_primersearch \
cB_printsextract \
cB_profit \
cB_prophecy \
cB_prophet \
cB_prosextract \
cB_pscan \
cB_psiphi \
cB_rebaseextract \
cB_recoder \
cB_redata \
cB_refseqget \
cB_remap \
cB_restover \
cB_restrict \
cB_revseq \
cB_seealso \
cB_seqcount \
cB_seqmatchall \
cB_seqret \
cB_seqretsetall \
cB_seqretsplit \
cB_seqxref \
cB_seqxrefget \
cB_servertell \
cB_showalign \
cB_showdb \
cB_showfeat \
cB_showorf \
cB_showpep \
cB_showseq \
cB_showserver \
cB_shuffleseq \
cB_sigcleave \
cB_silent \
cB_sirna \
cB_sixpack \
cB_sizeseq \
cB_skipredundant \
cB_skipseq \
cB_splitsource \
cB_splitter \
cB_stretcher \
cB_stssearch \
cB_supermatcher \
cB_syco \
cB_taxget \
cB_taxgetdown \
cB_taxgetrank \
cB_taxgetspecies \
cB_taxgetup \
cB_tcode \
cB_textget \
cB_textsearch \
cB_tfextract \
cB_tfm \
cB_tfscan \
cB_tmap \
cB_tranalign \
cB_transeq \
cB_trimest \
cB_trimseq \
cB_trimspace \
cB_twofeat \
cB_union \
cB_urlget \
cB_variationget \
cB_vectorstrip \
cB_water \
cB_whichdb \
cB_wobble \
cB_wordcount \
cB_wordfinder \
cB_wordmatch \
cB_wossdata \
cB_wossinput \
cB_wossname \
cB_wossoperation \
cB_wossoutput \
cB_wossparam \
cB_wosstopic \
cB_xmlget \
cB_xmltext \
cB_yank\
cB_hmmalign\
cB_hmmbuild\
cB_hmmemit\
cB_hmmfetch\
cB_hmmpress\
cB_hmmscan\
cB_hmmsearch\
cB_cmalign\
cB_cmbuild\
cB_cmcalibrate\
cB_cmfetch\
cB_cmpress\
cB_cmscan\
cB_cmsearch\
cB_bwa_idx\
cB_bwa_mem\
cB_dds\
cB_dps\
cB_ext\
cB_filter\
cB_gap2\
cB_nap\
cB_show

dist_pkgdata_DATA = \
emboss.default.template \
emboss.standard \
server.biomart \
server.dasregistry \
server.dbfetch server.dkfz \
server.ebeye \
server.entrez \
server.ensembl \
server.ensemblgenomes \
server.mrs \
server.mrs3 \
server.mrs4 \
server.srs \
server.wsdbfetch

# EMBOSS bin application SOURCES

cB_aaindexextract_SOURCES = aaindexextract.c
cB_abiview_SOURCES = abiview.c
cB_acdc_SOURCES = acdc.c
cB_acdgalaxy_SOURCES = acdgalaxy.c
cB_acdpretty_SOURCES = acdpretty.c
cB_acdtable_SOURCES = acdtable.c
cB_acdtrace_SOURCES = acdtrace.c
cB_acdvalid_SOURCES = acdvalid.c
cB_aligncopy_SOURCES = aligncopy.c
cB_aligncopypair_SOURCES = aligncopypair.c
cB_antigenic_SOURCES = antigenic.c
cB_assemblyget_SOURCES = assemblyget.c
cB_backtranambig_SOURCES = backtranambig.c
cB_backtranseq_SOURCES = backtranseq.c
cB_banana_SOURCES = banana.c
cB_biosed_SOURCES = biosed.c
cB_btwisted_SOURCES = btwisted.c
cB_cachedas_SOURCES = cachedas.c
cB_cachedbfetch_SOURCES = cachedbfetch.c
cB_cacheebeyesearch_SOURCES = cacheebeyesearch.c
cB_cacheensembl_SOURCES = cacheensembl.c
cB_cai_SOURCES = cai.c
cB_chaos_SOURCES = chaos.c
cB_charge_SOURCES = charge.c
cB_checktrans_SOURCES = checktrans.c
cB_chips_SOURCES = chips.c
cB_cirdna_SOURCES = cirdna.c
cB_codcmp_SOURCES = codcmp.c
cB_codcopy_SOURCES = codcopy.c
cB_coderet_SOURCES = coderet.c
cB_compseq_SOURCES = compseq.c
cB_cons_SOURCES = cons.c
cB_consambig_SOURCES = consambig.c
cB_cpgplot_SOURCES = cpgplot.c
cB_cpgreport_SOURCES = cpgreport.c
cB_cusp_SOURCES = cusp.c
cB_cutgextract_SOURCES = cutgextract.c
cB_cutseq_SOURCES = cutseq.c
cB_dan_SOURCES = dan.c
cB_dbiblast_SOURCES = dbiblast.c
cB_dbifasta_SOURCES = dbifasta.c
cB_dbiflat_SOURCES = dbiflat.c
cB_dbigcg_SOURCES = dbigcg.c
cB_dbtell_SOURCES = dbtell.c
cB_dbxcompress_SOURCES = dbxcompress.c
cB_dbxedam_SOURCES = dbxedam.c
cB_dbxfasta_SOURCES = dbxfasta.c
cB_dbxflat_SOURCES = dbxflat.c
cB_dbxgcg_SOURCES = dbxgcg.c
cB_dbxobo_SOURCES = dbxobo.c
cB_dbxreport_SOURCES = dbxreport.c
cB_dbxresource_SOURCES = dbxresource.c
cB_dbxstat_SOURCES = dbxstat.c
cB_dbxtax_SOURCES = dbxtax.c
cB_dbxuncompress_SOURCES = dbxuncompress.c
cB_degapseq_SOURCES = degapseq.c
cB_density_SOURCES = density.c
cB_descseq_SOURCES = descseq.c
cB_diffseq_SOURCES = diffseq.c
cB_distmat_SOURCES = distmat.c
cB_dotmatcher_SOURCES = dotmatcher.c
cB_dotpath_SOURCES = dotpath.c
cB_dottup_SOURCES = dottup.c
cB_dreg_SOURCES = dreg.c
cB_drfinddata_SOURCES = drfinddata.c
cB_drfindformat_SOURCES = drfindformat.c
cB_drfindid_SOURCES = drfindid.c
cB_drfindresource_SOURCES = drfindresource.c
cB_drget_SOURCES = drget.c
cB_drtext_SOURCES = drtext.c
cB_edamdef_SOURCES = edamdef.c
cB_edamhasinput_SOURCES = edamhasinput.c
cB_edamhasoutput_SOURCES = edamhasoutput.c
cB_edamisformat_SOURCES = edamisformat.c
cB_edamisid_SOURCES = edamisid.c
cB_edamname_SOURCES = edamname.c
cB_edialign_SOURCES = edialign.c
cB_einverted_SOURCES = einverted.c
cB_embossdata_SOURCES = embossdata.c
cB_embossversion_SOURCES = embossversion.c
cB_emma_SOURCES = emma.c
cB_emowse_SOURCES = emowse.c
cB_entret_SOURCES = entret.c
cB_epestfind_SOURCES = epestfind.c
cB_eprimer3_SOURCES = eprimer3.c
cB_eprimer32_SOURCES = eprimer32.c
cB_equicktandem_SOURCES = equicktandem.c
cB_est2genome_SOURCES = est2genome.c
cB_etandem_SOURCES = etandem.c
cB_extractalign_SOURCES = extractalign.c
cB_extractfeat_SOURCES = extractfeat.c
cB_extractseq_SOURCES = extractseq.c
cB_featcopy_SOURCES = featcopy.c
cB_featmerge_SOURCES = featmerge.c
cB_featreport_SOURCES = featreport.c
cB_feattext_SOURCES = feattext.c
cB_findkm_SOURCES = findkm.c
cB_freak_SOURCES = freak.c
cB_fuzznuc_SOURCES = fuzznuc.c
cB_fuzzpro_SOURCES = fuzzpro.c
cB_fuzztran_SOURCES = fuzztran.c
cB_garnier_SOURCES = garnier.c
cB_geecee_SOURCES = geecee.c
cB_getorf_SOURCES = getorf.c
cB_godef_SOURCES = godef.c
cB_goname_SOURCES = goname.c
cB_helixturnhelix_SOURCES = helixturnhelix.c
cB_hmoment_SOURCES = hmoment.c
cB_iep_SOURCES = iep.c
cB_infoalign_SOURCES = infoalign.c
cB_infoassembly_SOURCES = infoassembly.c
cB_infobase_SOURCES = infobase.c
cB_inforesidue_SOURCES = inforesidue.c
cB_infoseq_SOURCES = infoseq.c
cB_isochore_SOURCES = isochore.c
cB_jaspextract_SOURCES = jaspextract.c
cB_jaspscan_SOURCES = jaspscan.c
cB_jembossctl_SOURCES = jembossctl.c
cB_lindna_SOURCES = lindna.c
cB_listor_SOURCES = listor.c
cB_makenucseq_SORCES = makenucseq.c
cB_makeprotseq_SOURCES = makeprotseq.c
cB_marscan_SOURCES = marscan.c
cB_maskambignuc_SOURCES = maskambignuc.c
cB_maskambigprot_SOURCES = maskambigprot.c
cB_maskfeat_SOURCES = maskfeat.c
cB_maskseq_SOURCES = maskseq.c
cB_matcher_SOURCES = matcher.c
cB_megamerger_SOURCES = megamerger.c
cB_merger_SOURCES = merger.c
cB_msbar_SOURCES = msbar.c
cB_mwcontam_SOURCES = mwcontam.c
cB_mwfilter_SOURCES = mwfilter.c
cB_needle_SOURCES = needle.c
cB_needleall_SOURCES = needleall.c
cB_newcpgreport_SOURCES = newcpgreport.c
cB_newcpgseek_SOURCES = newcpgseek.c
cB_newseq_SOURCES = newseq.c
cB_nohtml_SOURCES = nohtml.c
cB_noreturn_SOURCES = noreturn.c
cB_nospace_SOURCES = nospace.c
cB_notab_SOURCES = notab.c
cB_notseq_SOURCES = notseq.c
cB_nthseq_SOURCES = nthseq.c
cB_nthseqset_SOURCES = nthseqset.c
cB_octanol_SOURCES = octanol.c
cB_oddcomp_SOURCES = oddcomp.c
cB_ontocount_SOURCES = ontocount.c
cB_ontoget_SOURCES = ontoget.c
cB_ontogetcommon_SOURCES = ontogetcommon.c
cB_ontogetdown_SOURCES = ontogetdown.c
cB_ontogetobsolete_SOURCES = ontogetobsolete.c
cB_ontogetroot_SOURCES = ontogetroot.c
cB_ontogetsibs_SOURCES = ontogetsibs.c
cB_ontogetup_SOURCES = ontogetup.c
cB_ontoisobsolete_SOURCES = ontoisobsolete.c
cB_ontotext_SOURCES = ontotext.c
cB_palindrome_SOURCES = palindrome.c
cB_pasteseq_SOURCES = pasteseq.c
cB_patmatdb_SOURCES = patmatdb.c
cB_patmatmotifs_SOURCES = patmatmotifs.c
cB_pepcoil_SOURCES = pepcoil.c
cB_pepdigest_SOURCES = pepdigest.c
cB_pepinfo_SOURCES = pepinfo.c
cB_pepnet_SOURCES = pepnet.c
cB_pepstats_SOURCES = pepstats.c
cB_pepwheel_SOURCES = pepwheel.c
cB_pepwindow_SOURCES = pepwindow.c
cB_pepwindowall_SOURCES = pepwindowall.c
cB_plotcon_SOURCES = plotcon.c
cB_plotorf_SOURCES = plotorf.c
cB_polydot_SOURCES = polydot.c
cB_preg_SOURCES = preg.c
cB_prettyplot_SOURCES = prettyplot.c
cB_prettyseq_SOURCES = prettyseq.c
cB_primersearch_SOURCES = primersearch.c
cB_printsextract_SOURCES = printsextract.c
cB_profit_SOURCES = profit.c
cB_prophecy_SOURCES = prophecy.c
cB_prophet_SOURCES = prophet.c
cB_prosextract_SOURCES = prosextract.c
cB_pscan_SOURCES = pscan.c
cB_psiphi_SOURCES = psiphi.c
cB_rebaseextract_SOURCES = rebaseextract.c
cB_recoder_SOURCES = recoder.c
cB_redata_SOURCES = redata.c
cB_refseqget_SOURCES = refseqget.c
cB_remap_SOURCES = remap.c
cB_restover_SOURCES = restover.c
cB_restrict_SOURCES = restrict.c
cB_revseq_SOURCES = revseq.c
cB_seealso_SOURCES = seealso.c
cB_seqcount_SOURCES = seqcount.c
cB_seqmatchall_SOURCES = seqmatchall.c
cB_seqret_SOURCES = seqret.c
cB_seqretsetall_SOURCES = seqretsetall.c
cB_seqretsplit_SOURCES = seqretsplit.c
cB_seqxref_SOURCES = seqxref.c
cB_seqxrefget_SOURCES = seqxrefget.c
cB_servertell_SOURCES = servertell.c
cB_showalign_SOURCES = showalign.c
cB_showdb_SOURCES = showdb.c
cB_showfeat_SOURCES = showfeat.c
cB_showorf_SOURCES = showorf.c
cB_showpep_SOURCES = showpep.c
cB_showseq_SOURCES = showseq.c
cB_showserver_SOURCES = showserver.c
cB_shuffleseq_SOURCES = shuffleseq.c
cB_sigcleave_SOURCES = sigcleave.c
cB_silent_SOURCES = silent.c
cB_sirna_SOURCES = sirna.c
cB_sixpack_SOURCES = sixpack.c
cB_sizeseq_SOURCES = sizeseq.c
cB_skipredundant_SOURCES = skipredundant.c
cB_skipseq_SOURCES = skipseq.c
cB_splitsource_SOURCES = splitsource.c
cB_splitter_SOURCES = splitter.c
cB_stretcher_SOURCES = stretcher.c
cB_stssearch_SOURCES = stssearch.c
cB_supermatcher_SOURCES = supermatcher.c
cB_syco_SOURCES = syco.c
cB_taxget_SOURCES = taxget.c
cB_taxgetdown_SOURCES = taxgetdown.c
cB_taxgetrank_SOURCES = taxgetrank.c
cB_taxgetspecies_SOURCES = taxgetspecies.c
cB_taxgetup_SOURCES = taxgetup.c
cB_tcode_SOURCES = tcode.c
cB_textget_SOURCE = textget.c
cB_textsearch_SOURCES = textsearch.c
cB_tfextract_SOURCES = tfextract.c
cB_tfm_SOURCES = tfm.c
cB_tfscan_SOURCES = tfscan.c
cB_tmap_SOURCES = tmap.c
cB_tranalign_SOURCES = tranalign.c
cB_transeq_SOURCES = transeq.c
cB_trimest_SOURCES = trimest.c
cB_trimseq_SOURCES = trimseq.c
cB_trimspace_SOURCES = trimspace.c
cB_twofeat_SOURCES = twofeat.c
cB_union_SOURCES = union.c
cB_urlget_SOURCES = urlget.c
cB_variationget_SOURCES = variationget.c
cB_vectorstrip_SOURCES = vectorstrip.c
cB_water_SOURCES = water.c
cB_whichdb_SOURCES = whichdb.c
cB_wobble_SOURCES = wobble.c
cB_wordcount_SOURCES = wordcount.c
cB_wordfinder_SOURCES = wordfinder.c
cB_wordmatch_SOURCES = wordmatch.c
cB_wossdata_SOURCES = wossdata.c
cB_wossinput_SOURCES = wossinput.c
cB_wossname_SOURCES = wossname.c
cB_wossoperation_SOURCES = wossoperation.c
cB_wossoutput_SOURCES = wossoutput.c
cB_wossparam_SOURCES = wossparam.c
cB_wosstopic_SOURCES = wosstopic.c
cB_xmlget_SOURCES = xmlget.c
cB_xmltext_SOURCES = xmltext.c
cB_yank_SOURCES = yank.c
cB_hmmalign_SOURCES = hmmalign.c
cB_hmmbuild_SOURCES = hmmbuild.c
cB_hmmemit_SOURCES = hmmemit.c
cB_hmmfetch_SOURCES = hmmfetch.c
cB_hmmpress_SOURCES = hmmpress.c
cB_hmmscan_SOURCES = hmmscan.c
cB_hmmsearch_SOURCES = hmmsearch.c
cB_cmalign_SOURCES = cmalign.c
cB_cmbuild_SOURCES = cmbuild.c
cB_cmcalibrate_SOURCES = cmcalibrate.c
cB_cmfetch_SOURCES = cmfetch.c
cB_cmpress_SOURCES = cmpress.c
cB_cmscan_SOURCES = cmscan.c
cB_cmsearch_SOURCES = cmsearch.c
cB_bwa_idx_SOURCES = bwa-idx.c
cB_bwa_mem_SOURCES = bwa-mem.c
cB_dds_SOURCES = dds.c
cB_dps_SOURCES = dps.c
cB_ext_SOURCES = ext.c
cB_filter_SOURCES = filter.c
cB_gap2_SOURCES = gap2.c
cB_nap_SOURCES = nap.c
cB_show_SOURCES = show.c
cB_embossupdate_SOURCES = embossupdate.c 
cB_makenucseq_SOURCES = makenucseq.c
cB_textget_SOURCES = textget.c

# EMBOSS check application SOURCES

cB_acdlog_SOURCES = acdlog.c
cB_acdrelations_SOURCES = acdrelations.c
cB_acdxsd_SOURCES = acdxsd.c
cB_ajtest_SOURCES = ajtest.c
cB_complex_SOURCES = complex.c
cB_corbatest_SOURCES = corbatest.c
cB_dastest_SOURCES = dastest.c
cB_dbfetchtest_SOURCES = dbfetchtest.c
cB_domtesta_SOURCES = domtesta.c
cB_domtestb_SOURCES = domtestb.c
cB_domtestc_SOURCES = domtestc.c
cB_domtestd_SOURCES = domtestd.c
cB_edamclean_SOURCES = edamclean.c
cB_ensembltest_SOURCES = ensembltest.c
cB_entrails_SOURCES = entrails.c
cB_entrailsbook_SOURCES = entrailsbook.c
cB_entrailshtml_SOURCES = entrailshtml.c
cB_entrailswiki_SOURCES = entrailswiki.c
cB_histogramtest_SOURCES = histogramtest.c
cB_intconv_SOURCES = intconv.c
cB_martattributes_SOURCES = martattributes.c
cB_martdatasets_SOURCES = martdatasets.c
cB_martfilters_SOURCES = martfilters.c
cB_martquery_SOURCES = martquery.c
cB_martregistry_SOURCES = martregistry.c
cB_martseqs_SOURCES = martseqs.c
cB_newcoils_SOURCES = newcoils.c
cB_origsplitter_SOURCES = origsplitter.c
cB_origunion_SOURCES = origunion.c
cB_patmattest_SOURCES = patmattest.c
cB_prima_SOURCES = prima.c
cB_primers_SOURCES = primers.c
cB_seqinfo_SOURCES = seqinfo.c
cB_seqretall_SOURCES = seqretall.c
cB_seqretallfeat_SOURCES = seqretallfeat.c
cB_seqretset_SOURCES = seqretset.c
cB_seqretsingle_SOURCES = seqretsingle.c
cB_seqrettype_SOURCES = seqrettype.c
cB_sqltest_SOURCES = sqltest.c
cB_testplot_SOURCES = testplot.c
cB_treetypedisplay_SOURCES = treetypedisplay.c
