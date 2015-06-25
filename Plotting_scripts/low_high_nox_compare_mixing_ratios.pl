#! /usr/bin/env perl
# Compare mixing ratio time series of species listed in ARGV over different NOx conditions, non-tagged species only
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use Statistics::R;

my @runs = qw( NOx_limited NOx_saturated Tuned );
my @mechanisms = qw( MOZART-4 );
my $base = "/work/users/jco/Variable_Conditions";
my %data;

my $mecca = MECCA->new("$base/Tuned/MOZART-4/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    foreach my $run (@runs) {
        my $dir = "$base/$run/$mechanism";
        my $mecca = MECCA->new("$dir/boxmodel");
        foreach my $spc (qw(CO OH HO2 O3 NO NO2)) {
            my $mixing_ratio = $mecca->tracer($spc);
            $data{$mechanism}{$run}{$spc} = $mixing_ratio(1:$ntime-2);
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

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` d = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $run (sort keys %{$data{$mechanism}}) {
        $R->set('run', $run);
        $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
        foreach my $spc (sort keys %{$data{$mechanism}{$run}}) {
            $R->set('species', $spc);
            $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$run}{$spc}->dog ]);
            $R->run(q` pre[species] = mixing.ratio `);
        }
        $R->run(q` pre = gather(pre, Species, Mixing.Ratio, -Time, -Run) `,
                q` d = rbind(d, pre) `,
        );
    }
}
$R->run(q` d$Species = factor(d$Species, levels = c("O3", "OH", "HO2", "CO", "NO", "NO2")) `);
#my $p = $R->run(q` print(d) `);
#print $p, "\n";

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, colour = Run, order = Run, group = Run)) `,
        q` p = p + geom_line(size = 1) `,
        q` p = p + facet_wrap( ~ Species, scales = "free_y") `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` p = p + xlab("Time (days)") + ylab("Mixing Ratio") `, 
        q` p = p + ggtitle("Mixing Ratio Comparisons in Different Atmospheric Regimes") `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(strip.text = element_text(size = 15, face = "bold")) `,
        q` p = p + theme(plot.title = element_text(size = 18, face = "bold")) `,
        q` p = p + theme(axis.title = element_text(size = 14, face = "bold")) `,
        q` p = p + theme(axis.text = element_text(size = 13)) `,
        q` p = p + theme(legend.text = element_text(size = 13)) `,
        q` p = p + theme(panel.margin = unit(5, "mm")) `,
        q` p = p + theme(legend.position = "top") `,
);

$R->run(q` CairoPDF(file = "low_high_nox_non-tagged_mixing_ratios_comparison.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
