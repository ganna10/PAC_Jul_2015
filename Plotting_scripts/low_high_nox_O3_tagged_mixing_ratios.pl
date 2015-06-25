#! /usr/bin/env perl
# Time series of all tagged O3 mixing ratios stacked to give total O3, facetted by different regime
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $species = "O3";
my @runs = qw( NOx_limited NOx_saturated Tuned );
my @mechanisms = qw( MOZART-4 );
my $base = "/work/users/jco/Variable_Conditions";
my %data;

my $mecca = MECCA->new("$base/Tuned/MOZART-4/boxmodel");
my $ntime = $mecca->time->nelem;
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;
$time = $time(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    foreach my $run (@runs) {
        my $dir = "$base/$run/$mechanism";
        my $mecca = MECCA->new("$dir/boxmodel");
        my @tagged_species = get_tagged_species($dir, $species);
        foreach my $spc (@tagged_species) {
            my ($base, $tag) = split /_X_/, $spc;
            my $mixing_ratio = $mecca->tracer($spc);
            $mixing_ratio = $mixing_ratio(1:$ntime-2);
            $data{$mechanism}{$run}{$tag} = $mixing_ratio;
        }
    }
}

foreach my $mechanism (keys %data) { #allocate MOZART species to Functional Groups
    foreach my $run (keys %{$data{$mechanism}}) {
        foreach my $spc (keys %{$data{$mechanism}{$run}}) {
            if ($spc eq "C2H6" or $spc eq "C3H8" or $spc eq "BIGALK") {
                $data{$mechanism}{$run}{"Alkanes"} += $data{$mechanism}{$run}{$spc};
                delete $data{$mechanism}{$run}{$spc};
            } elsif ($spc eq "C2H4" or $spc eq "C3H6" or $spc eq "BIGENE") {
                $data{$mechanism}{$run}{"Alkenes"} += $data{$mechanism}{$run}{$spc};
                delete $data{$mechanism}{$run}{$spc};
            } elsif ($spc eq "TOLUENE") {
                $data{$mechanism}{$run}{"Aromatics"} += $data{$mechanism}{$run}{$spc};
                delete $data{$mechanism}{$run}{$spc};
            } elsif ($spc eq "ISOP") {
                $data{$mechanism}{$run}{"Isoprene"} += $data{$mechanism}{$run}{$spc};
                delete $data{$mechanism}{$run}{$spc};
            }
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);

$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` d = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $run (sort keys %{$data{$mechanism}}) {
        $R->set('run', $run);
        $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
        foreach my $voc (sort keys %{$data{$mechanism}{$run}}) {
            $R->set('voc', $voc);
            $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$run}{$voc}->dog ]);
            $R->run(q` pre[voc] = mixing.ratio * 1e9 `);
        }
        $R->run(q` pre = gather(pre, VOC, Mixing.Ratio, -Time, -Run) `,
                q` d = rbind(d, pre) `,
        );
    }
}
$R->run(q` d$VOC = factor(d$VOC, levels = c("INI", "XTR", "CO", "CH4", "Alkanes", "Alkenes", "Isoprene", "Aromatics")) `,
        q` my.colours = c("INI" = "#6c254f", "XTR" = "#f9c500", "CO" = "#0e5c28", "CH4" = "#2b9eb3", "Alkanes" = "#ef6638", "Alkenes" = "#86c650", "Isoprene" = "#8c1531", "Aromatics" = "#0352cb") `,
);

$R->run(q` d$Run = factor(d$Run, levels = c("NOx_limited", "Tuned", "NOx_saturated")) `,
        q` run.colours = c("NOx_limited" = "#6c254f", "Tuned" = "#2b9eb3", "NOx_saturated" = "#ef6638") `,
);

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, fill = VOC, order = VOC)) `,
        q` p = p + geom_area(position = "stack") `,
        q` p = p + geom_area(position = "stack", colour = "black", show_guide = FALSE) `,
        q` p = p + facet_wrap( ~ Run) `,
        q` p = p + theme_tufte() `,
        q` p = p + ggtitle("Allocated Ozone Mixing Ratios") `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` p = p + scale_y_continuous(expand = c(0, 1e-9)) `,
        q` p = p + ylab("Mixing Ratio (ppbv)") `,
        q` p = p + theme(strip.text = element_text(size = 15, face = "bold")) `,
        q` p = p + theme(plot.title = element_text(size = 18, face = "bold")) `,
        q` p = p + theme(axis.title = element_text(size = 14, face = "bold")) `,
        q` p = p + theme(axis.text = element_text(size = 13)) `,
        q` p = p + theme(legend.text = element_text(size = 13)) `,
        q` p = p + theme(panel.margin = unit(5, "mm")) `,
        q` p = p + xlab("Time (Days)") `,
        q` p = p + theme(legend.position = "top") `,
        q` p = p + scale_fill_manual(values = my.colours, limits = rev(levels(d$VOC))) `,
);

$R->run(q` CairoPDF(file = "low_high_nox_O3_mixing_ratio_components_facet_run.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, colour = Run, order = Run)) `,
        q` p = p + geom_line(size = 1) `,
        q` p = p + facet_wrap( ~ VOC, scales = "free_y", nrow = 2) `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        #q` p = p + scale_y_continuous(expand = c(0, 1e-9)) `,
        q` p = p + theme(strip.text = element_text(face = "bold")) `,
        q` p = p + ylab("Mixing Ratio (ppbv)") `,
        q` p = p + xlab("Time (Days)") `,
        q` p = p + theme(legend.position = "top") `,
        q` p = p + scale_fill_manual(values = run.colours, limits = rev(levels(d$Run))) `,
);

$R->run(q` CairoPDF(file = "low_high_nox_O3_mixing_ratio_components_facet_group.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_tagged_species {
    my ($dir, $species) = @_;
    my $spc_file = "$dir/gas.spc";
    my @tagged;
    open my $in, '<:encoding(utf-8)', $spc_file or die $!;
    while (<$in>) {
        next unless ($_ =~ /^${species}_/);
        my ($spc, $rest) = split / = /, $_;
        push @tagged, $spc;
    }
    close $in;
    return @tagged;
}

sub get_voc_name {
    my ($spc) = @_;
    my $VOC;
    if ($spc eq "C2H6") {
        $VOC = "Ethane";
    } elsif ($spc eq "C3H8") {
        $VOC = "Propane";
    } elsif ($spc eq "NC4H10") {
        $VOC = "Butane";
    } elsif ($spc eq "IC4H10") {
        $VOC = "2-Methylpropane";
    } elsif ($spc eq "NC5H12") {
        $VOC = "Pentane";
    } elsif ($spc eq "IC5H12") {
        $VOC = "2-Methylbutane";
    } elsif ($spc eq "NC6H14") {
        $VOC = "Hexane";
    } elsif ($spc eq "NC7H16") {
        $VOC = "Heptane";
    } elsif ($spc eq "NC8H18") {
        $VOC = "Octane";
    } elsif ($spc eq "C2H4") {
        $VOC = "Ethene";
    } elsif ($spc eq "C3H6") {
        $VOC = "Propene";
    } elsif ($spc eq "BUT1ENE") {
        $VOC = "Butene";
    } elsif ($spc eq "MEPROPENE") {
        $VOC = "2-Methylpropene";
    } elsif ($spc eq "ISOP") {
        $VOC = "Isoprene";
    } elsif ($spc eq "BENZENE") {
        $VOC = "Benzene";
    } elsif ($spc eq "TOLUENE" or $spc eq "TOLUENE_M") {
        $VOC = "Toluene";
    } elsif ($spc eq "MXYL") {
        $VOC = "m-Xylene";
    } elsif ($spc eq "OXYL") {
        $VOC = "o-Xylene";
    } elsif ($spc eq "PXYL") {
        $VOC = "p-Xylene";
    } elsif ($spc eq "EBENZ") {
        $VOC = "Ethylbenzene";
    } elsif ($spc eq "CH4") {
        $VOC = "Methane";
    } elsif ($spc eq "XTR" or $spc eq "INI" or "CO") {
        $VOC = $spc;
    } else {
        print "No mapping for $spc\n";
    }
    return $VOC;
}
