/* HTx follow-up table: Excel results to RTF via PROC REPORT.
   Creates a separate TableHeader file and a table with a footer.
   Expected Excel columns: row_order, row_id, section, row_type,
   label, indent, footnote, p0_7, ge7, entire. */

/* Paths and report settings */
%let indir = D:\RNC\SCAARPCAV\Output;
%let outdir = D:\RNC\SCAARPCAV\Output;
%let infile = &indir.\htx_followup_table_for_sas.xlsx;
%let FILE_NUMB = 01_01;
%let rdat8 = %sysfunc(today(),yymmddn8.);
%let font = Arial;
%let table_title = Table 1. Patients and follow-up, coronary angiography and PCI;

options orientation=landscape papersize=A4 nodate nonumber
        leftmargin=0.60in rightmargin=0.60in
        topmargin=0.65in bottommargin=0.65in;
ods escapechar='^';

/* Import the Excel table */
proc import datafile="&infile."
    out=work.htx_import dbms=xlsx replace;
    getnames=yes;
run;

/* Prepare text fields and numeric row order */
data work.htx_prepared;
    set work.htx_import;
    length label_print $600 rtype $30
           p0_7_print ge7_print entire_print $200
           order_num 8;

    order_num = input(
        tranwrd(strip(vvalue(row_order)), ',', '.'), best32.
    );
    rtype = lowcase(strip(row_type));
    label_print = strip(label);

    /* Replace UTF-8 comparison signs with ASCII, if still present.
       If Excel import already converted them to ?, replace in R
       before exporting; the original symbol cannot be recovered here. */
    label_print = tranwrd(label_print, 'E289A5'x, '>=');
    label_print = tranwrd(label_print, 'E289A4'x, '<=');

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

/* Carry section titles forward to data rows for PROC REPORT */
data work.htx_report;
    set work.htx_prepared;
    length section_title $300;
    retain section_title;
    if rtype = 'section' then do;
        section_title = strip(label_print);
        delete;
    end;
run;

/* Separate title RTF, following the usual TableHeader convention */
options orientation=landscape;
title; footnote;
ods rtf file="&outdir.\TableHeader_&FILE_NUMB._&rdat8..rtf"
    style=styles.journal bodytitle;
ods rtf text="^S={font=('Arial',11pt,bold)}&table_title.";
ods rtf close;

/* Main RTF: title is deliberately omitted because it is in TableHeader */
ods rtf file="&outdir.\Table_&FILE_NUMB._&rdat8..rtf"
    style=styles.statistical bodytitle startpage=no;

proc report data=work.htx_report nowd missing
    style(report)=[frame=box rules=all
                   bordercolor=cxB8B8B8 borderwidth=0.5pt
                   cellpadding=3pt cellspacing=0]
    style(header)=[background=cxC4C4C4 foreground=black
                   fontfamily="&font" fontsize=9pt
                   fontweight=bold just=center vjust=middle]
    style(column)=[background=white foreground=black
                   fontfamily="&font" fontsize=9pt vjust=middle];

    column section_title order_num rtype indent
           label_print p0_7_print ge7_print entire_print;

    define section_title / order order=data noprint;
    define order_num / order noprint;
    define rtype / display noprint;
    define indent / display noprint;

    /* A4 landscape: one wide label column, three equal result columns */
    define label_print / display 'Characteristic'
        style(column)=[cellwidth=4.10in just=left]
        style(header)=[just=center];
    define p0_7_print / display 'Follow-up 0-<7 years'
        style(column)=[cellwidth=2.05in just=center];
    define ge7_print / display 'Follow-up >=7 years'
        style(column)=[cellwidth=2.05in just=center];
    define entire_print / display 'Entire follow-up'
        style(column)=[cellwidth=2.05in just=center];

    /* Section titles span all four displayed columns */
    compute before section_title /
        style=[background=cxDADADA foreground=black
               fontfamily="&font" fontsize=9pt
               fontweight=bold just=center];
        line section_title $300.;
    endcomp;

    /* Indent levels within the RTF cell; do not rely on spaces */
    compute label_print;
        if indent > 0 then
            call define(_col_, 'style/merge',
                cats('style=[indent=',
                     strip(put(0.22 * indent, 5.2)), 'in]'));
        if rtype = 'subheading' then
            call define(_col_, 'style/merge',
                'style=[fontweight=bold]');
    endcomp;

    /* Explanatory footer below the final table row */
    compute after / style=[background=cxECEDEC just=left
                           fontfamily="&font" fontsize=8pt];
        line 'Data are presented as n (%), median (IQR), or rates, as indicated.';
        line 'Follow-up periods: 0 to <7 years, >=7 years, and entire follow-up.';
    endcomp;
run;

ods rtf close;
title; footnote;

%put NOTE: Title RTF: &outdir.\TableHeader_&FILE_NUMB._&rdat8..rtf;
%put NOTE: Table RTF: &outdir.\Table_&FILE_NUMB._&rdat8..rtf;
