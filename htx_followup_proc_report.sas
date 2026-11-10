/* HTx follow-up table: Excel -> SAS ODS RTF using PROC REPORT.
   Excel columns: row_order row_id section row_type label indent footnote p0_7 ge7 entire.
   Section rows are converted to merged grey bands; only data/subheading rows print. */

/* File paths */
%let indir = D:\RNC\SCAARPCAV\Output;
%let outdir = D:\RNC\SCAARPCAV\Output;
%let infile = &indir.\htx_followup_table_for_sas.xlsx;
%let rdat8 = %sysfunc(today(),yymmddn8.);
%let outfile = &outdir.\Table_1_&rdat8..rtf;
%let font = Arial;  /* Switch to Aptos if installed on the SAS host */

options orientation=portrait papersize=A4 nodate nonumber
        leftmargin=0.60in rightmargin=0.60in
        topmargin=0.65in bottommargin=0.65in;
ods escapechar='^';

/* Import Excel */
proc import datafile="&infile."
    out=work.htx_import dbms=xlsx replace;
    getnames=yes;
run;

/* Convert mixed Excel types to printable strings */
data work.htx_prepared;
    set work.htx_import;
    length label_print $600 rtype $30
           p0_7_print ge7_print entire_print $200
           order_num 8;
    order_num = input(tranwrd(strip(vvalue(row_order)), ',', '.'), best32.);
    rtype = lowcase(strip(row_type));
    label_print = strip(label);
    if indent > 0 then label_print = cats('^{nbspace 3}', label_print);
    if lowcase(strip(vvalue(footnote))) = 'c' then
        label_print = cats(label_print, '^{super c}');

    p0_7_print = strip(vvalue(p0_7));
    ge7_print = strip(vvalue(ge7));
    entire_print = strip(vvalue(entire));
    if missing(p0_7) then p0_7_print = ' ';
    if missing(ge7) then ge7_print = ' ';
    if missing(entire) then entire_print = ' ';
run;

proc sort data=work.htx_prepared;
    by order_num;
run;

/* Read the section titles from the special rows in Excel.
   Carry the most recent title forward to ordinary rows. */
data work.htx_report;
    set work.htx_prepared;
    length section_title $300;
    retain section_title;
    if rtype = 'section' then do;
        section_title = strip(label);
        delete;
    end;
run;

/* Main table, entirely with PROC REPORT */
ods rtf file="&outfile." style=styles.statistical startpage=no;

title1 j=l font="&font" height=10pt
    'Table 1. Patients and follow-up, coronary angiography and PCI';

/* Add your final study-specific denominator definitions here, if needed.
footnote1 j=l font="&font" height=8pt
    'Percentages are calculated using the stated denominators.';
*/

proc report data=work.htx_report nowd missing
    style(report)=[frame=box rules=all bordercolor=cxB8B8B8
                   borderwidth=0.5pt cellpadding=3pt cellspacing=0]
    style(header)=[background=cxC4C4C4 foreground=black
                   fontfamily="&font" fontsize=9pt fontweight=bold
                   just=center vjust=middle]
    style(column)=[background=white foreground=black
                   fontfamily="&font" fontsize=9pt vjust=middle];

    column section_title order_num rtype
           label_print p0_7_print ge7_print entire_print;

    /* Hidden ordering and formatting metadata */
    define section_title / order order=data noprint;
    define order_num / order noprint;
    define rtype / display noprint;

    /* Four visible columns */
    define label_print / display 'Characteristic'
        style(column)=[cellwidth=3.0in just=left]
        style(header)=[just=center];
    define p0_7_print / display 'Follow-up 0-<7 years'
        style(column)=[cellwidth=1.55in just=center];
    define ge7_print / display 'Follow-up >=7 years'
        style(column)=[cellwidth=1.10in just=center];
    define entire_print / display 'Entire follow-up'
        style(column)=[cellwidth=1.10in just=center];

    /* Section headings become merged, grey rows across the report */
    compute before section_title /
        style=[background=cxDADADA foreground=black
               fontfamily="&font" fontsize=9pt fontweight=bold
               just=center bordertopcolor=cxB8B8B8
               borderbottomcolor=cxB8B8B8];
        line section_title $300.;
    endcomp;

    /* Emphasize the intermediate subheading, but not as a new section */
    compute label_print;
        if rtype = 'subheading' then
            call define('label_print', 'style', 'style=[fontweight=bold]');
    endcomp;
run;

ods rtf close;
title; footnote;
%put NOTE: Created &outfile.;
