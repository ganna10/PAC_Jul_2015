#! /usr/bin/env perl
# Vertical mixing setup mixing ratios
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use MECCA;
use Statistics::R;

my $dir = "/work/users/jco/Meteorology_Tests/Vertical_Mixing";
my $mecca = MECCA->new("$dir/boxmodel");
my @species = qw(O3 OH HO2 NO NO2 CO);
my %data;

my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;

foreach my $spc (@species) {
    my $mixing_ratio = $mecca->tracer($spc);
    $data{$spc} = $mixing_ratio;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` d = data.frame(Time) `);
foreach my $spc (sort keys %data) {
    $R->set('species', $spc);
    $R->set('mixing.ratio', [ map { $_ } $data{$spc}->dog ]);
    $R->run(q` d[species] = mixing.ratio `);
}
$R->run(q` d = gather(d, Species, Mixing.Ratio, -Time) `);
$R->run(q` d$Species = factor(d$Species, levels = c("O3", "OH", "HO2", "CO", "NO", "NO2")) `);

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio)) `,
        q` p = p + geom_line(size = 1) `,
        q` p = p + facet_wrap( ~ Species, scales = "free_y") `,
        q` p = p + theme_tufte() `,
        q` p = p + ggtitle("Mixing Ratios in Vertical Mixing set-up") `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` p = p + ylab("Mixing Ratio") `,
        q` p = p + xlab("Time (Days)") `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(strip.text = element_text(size = 15, face = "bold")) `,
        q` p = p + theme(plot.title = element_text(size = 18, face = "bold")) `,
        q` p = p + theme(axis.title = element_text(size = 14, face = "bold")) `,
        q` p = p + theme(axis.text = element_text(size = 13)) `,
        q` p = p + theme(legend.text = element_text(size = 13)) `,
        q` p = p + theme(panel.margin = unit(5, "mm")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + theme(legend.position = "top") `,
);

$R->run(q` CairoPDF(file = "vertical_mixing_mixing_ratios.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
