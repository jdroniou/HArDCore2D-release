#!/usr/bin/perl

use Template;

my @errors = (
    "ErrUL2",
    "ErrURotRot",
    "ErrPL2",
    "ErrPGrad"
    );

my $vars = {
    errors => \@errors
};

my $template = Template->new();
$template->process("runseries.tmpl", $vars, "runseries.sh") or die $template->error();
system("chmod +x runseries.sh")
