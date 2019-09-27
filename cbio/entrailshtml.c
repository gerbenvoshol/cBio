/* @source entrails application
**
** Show EMBOSS internals
**
** @author Copyright (C) Peter Rice
** @@
**
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
******************************************************************************/
#include "cbio.h"




/* @prog entrailshtml *********************************************************
**
** Reports internal data structures in HTML
**
******************************************************************************/

int main(int argc, char **argv)
{

	AjPFile outf;
	AjBool doall = ajFalse;
	AjPStr *tables = NULL;
	ajuint i;

	embInit("entrailshtml", argc, argv);

	outf = ajAcdGetOutfile("outfile");
	tables = ajAcdGetList("tables");

	for (i = 0; tables[i]; i++) {
		if (ajStrMatchC(tables[i], "all")) {
			doall = ajTrue;
		}

		if (doall || ajStrMatchC(tables[i], "inseq")) {
			ajSeqPrinthtmlInFormat(outf);
		}

		if (doall || ajStrMatchC(tables[i], "outseq")) {
			ajSeqoutPrinthtmlFormat(outf);
		}

		if (doall || ajStrMatchC(tables[i], "infeat")) {
			ajFeatinPrinthtmlFormat(outf);
		}

		if (doall || ajStrMatchC(tables[i], "outfeat")) {
			ajFeatoutPrinthtmlFormat(outf);
		}

		if (doall || ajStrMatchC(tables[i], "report")) {
			ajReportPrinthtmlFormat(outf);
		}

		if (doall || ajStrMatchC(tables[i], "align")) {
			ajAlignPrinthtmlFormat(outf);
		}

		if (doall || ajStrMatchC(tables[i], "assembly")) {
			ajAsseminprintHtml(outf);
			ajAssemoutprintHtml(outf);
		}

		if (doall || ajStrMatchC(tables[i], "resource")) {
			ajResourceinprintHtml(outf);
			ajResourceoutprintHtml(outf);
		}

		if (doall || ajStrMatchC(tables[i], "taxonomy")) {
			ajTaxinprintHtml(outf);
			ajTaxoutprintHtml(outf);
		}
	}

	ajFileClose(&outf);

	embExit();

	return 0;
}
