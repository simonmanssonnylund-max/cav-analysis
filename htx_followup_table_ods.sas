/* HTx follow-up table: Excel results -> formatted Word/RTF table
   Expected columns from R: row_order, row_id, section, row_type,
   label, indent, footnote, p0_7, ge7, entire.
   Uses DATA-step ODS (no PROC REPORT) so section rows span all 4 columns. */

/* 1. Paths and page settings */
%let infile  = D:\RNC\SCAARPCAV\Output\htx_followup_table_for_sas.xlsx;
%let outfile = D:\RNC\SCAARPCAV\Output\htx_followup_table.rtf;

options papersize=A4 orientation=portrait nodate nonumber
        leftmargin=0.60in rightmargin=0.60in
        topmargin=0.65in bottommargin=0.65in;
ods escapechar='^';

/* 2. Import the Excel sheet; by default PROC IMPORT reads the first sheet */
proc import datafile="&infile."
    out=work.htx_import
    dbms=xlsx replace;
    getnames=yes;
run;

/* Inspect imported column names/types if needed */
proc contents data=work.htx_import varnum; run;

/* 3. Prepare presentation strings and preserve row ordering */
data work.htx_table;
    length label_print $600;
    set work.htx_import;

    label_print = strip(label);
    if indent = 1 then label_print = cats('^{nbspace 3}', strip(label_print));
    if lowcase(strip(footnote)) = 'c' then
        label_print = cats(strip(label_print), '^{super c}');

    /* Missing formatted results should display as blank cells */
    if missing(p0_7) then p0_7 = ' ';
    if missing(ge7) then ge7 = ' ';
    if missing(entire) then entire = ' ';
run;

proc sort data=work.htx_table;
    by row_order;
run;

/* 4. Produce an editable RTF table using ODS */
ods rtf file="&outfile." style=styles.rtf startpage=no;

title1 j=l font='Aptos' height=10pt
    'Table 1. Patients, angiography and PCI during follow-up';

/* Add the footnote definitions here when finalized, for example:
footnote1 j=l font='Aptos' height=8pt '^{unicode 2020} [Diagnosis denominator definition].';
footnote2 j=l font='Aptos' height=8pt '^{super c} [PCI treatment definition].';
*/

data _null_;
    set work.htx_table end=last;
    length header_style section_style label_style number_style $500;
    length hdr1 hdr2 hdr3 hdr4 num1 num2 num3 $500;
    retain header_style section_style label_style number_style
           hdr1 hdr2 hdr3 hdr4 num1 num2 num3;

    if _n_ = 1 then do;
        declare odsout t();

        header_style = 'fontfamily=Aptos fontsize=9pt fontweight=bold ' ||
                       'color=black backgroundcolor=cxC4C4C4 ' ||
                       'bordercolor=cxB8B8B8 borderwidth=0.5pt';
        section_style = 'fontfamily=Aptos fontsize=9pt fontweight=bold ' ||
                        'color=black backgroundcolor=cxDADADA ' ||
                        'bordercolor=cxB8B8B8 borderwidth=0.5pt';
        label_style = 'fontfamily=Aptos fontsize=9pt color=black ' ||
                      'backgroundcolor=white bordercolor=cxB8B8B8 borderwidth=0.5pt ' ||
                      'cellwidth=2.85in';
        number_style = 'fontfamily=Aptos fontsize=9pt color=black ' ||
                       'backgroundcolor=white bordercolor=cxB8B8B8 borderwidth=0.5pt';

        hdr1 = cats(header_style, ' cellwidth=2.85in');
        hdr2 = cats(header_style, ' cellwidth=1.62in');
        hdr3 = cats(header_style, ' cellwidth=1.20in');
        hdr4 = hdr3;
        num1 = cats(number_style, ' cellwidth=1.62in');
        num2 = cats(number_style, ' cellwidth=1.20in');
        num3 = num2;

        t.table_start();
        t.head_start();
        t.row_start();
        t.format_cell(data:'Characteristic', just:'C',
            style_attr:hdr1);
        t.format_cell(data:'Follow-up 0-<7 years', just:'C',
            style_attr:hdr2);
        t.format_cell(data:'Follow-up >=7 years', just:'C',
            style_attr:hdr3);
        t.format_cell(data:'Entire follow-up', just:'C',
            style_attr:hdr4);
        t.row_end();
        t.head_end();
        t.body_start();
    end;

    t.row_start();
    if lowcase(strip(row_type)) = 'section' then do;
        /* A true merged row across the entire table */
        t.format_cell(data:label_print, column_span:4, just:'C',
            style_attr:section_style);
    end;
    else do;
        t.format_cell(data:label_print, just:'L', style_attr:label_style);
        t.format_cell(data:p0_7, just:'C',
            style_attr:num1);
        t.format_cell(data:ge7, just:'C',
            style_attr:num2);
        t.format_cell(data:entire, just:'C',
            style_attr:num3);
    end;
    t.row_end();

    if last then do;
        t.body_end();
        t.table_end();
    end;
run;

ods rtf close;
title; footnote;

%put NOTE: Saved the formatted table to &outfile.;
