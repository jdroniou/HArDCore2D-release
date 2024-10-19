#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV qw( csv );

open(my $fh, '>', 'table.tex');

my $directory = $ARGV[0] or die "A directory name must be provided\n";

# Table header
print $fh <<EOF;
\\documentclass{standalone}
\\usepackage{amsmath,amsfonts}
\\usepackage{booktabs}

\\begin{document}
\\begin{tabular}{cccccccccccc}
  \\toprule
  \$h\$
  & Spaces size
  & System size
  & \$L^2\$-error on \$\\boldsymbol{u}\$ & OCV
  & Energy error on \$\\boldsymbol{u}\$ & OCV
  & \$L^2\$-error on \$p\$ & OCV
  & \$H^1\$-error on \$p\$ & OCV
  \\\\
  \\midrule
EOF

# Table body
foreach my $k (0..2) {
    print $fh "\\multicolumn{11}{c}{\$ k = $k \$} \\\\ \n";    
    print $fh "\\midrule\n";
    
    my $input_file = "$directory/cart_k$k/data_rates.dat";
    my $data = csv(in => $input_file, headers => "auto");        
    foreach (@$data) {
        print $fh "$_->{'MeshSize'} & $_->{'DimSpaces'} & $_->{'DimLinSys'} & $_->{'ErrUL2'} & $_->{'ErrUL2_rate'} & $_->{'ErrURotRot'} & $_->{'ErrURotRot_rate'} & $_->{'ErrPL2'} & $_->{'ErrPL2_rate'} & $_->{'ErrPGrad'} & $_->{'ErrPGrad_rate'} \\\\ \n";
    }
    print $fh "\\midrule\n";
}

# Table footer
print $fh <<EOF;
  \\bottomrule
\\end{tabular}
\\end{document}
EOF
