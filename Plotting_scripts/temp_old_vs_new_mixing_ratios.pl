#! /usr/bin/env perl
# Compare mixing ratios between previous tagged runs and new tagged runs: 'real' chemistry only
# Version 0: Jane Coates 29/6/2015

use strict;
use diagnostics;
use MECCA;
use Statistics::R;

my $dir = "/work/users/jco/Meteorology_Tests/Temp_295";
my $old_tag_dir = "/work/users/jco/Meteorology_Tests/Temp_293";
my $old_mecca = MECCA->new("$old_tag_dir/boxmodel");
my $new_mecca = MECCA->new("$dir/boxmodel");
my %data;

my $time = $new_mecca->time;
$time -= $time->at(0);
$time /= 86400;

foreach my $spc (qw( OH O3 HO2 NO NO2 CO )) {
    $data{"295 K"}{$spc} += $new_mecca->tracer($spc);
    $data{"293 K"}{$spc} += $old_mecca->tracer($spc);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);
$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` d = data.frame() `);

foreach my $run (sort keys %data) {
    $R->set('run', $run);
    $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
    foreach my $spc (sort keys %{$data{$run}}) {
        $R->set('spc', $spc);
        $R->set('mixing.ratio', [ map { $_ } $data{$run}{$spc}->dog ]);
        $R->run(q` pre[spc] = mixing.ratio `);
    }
    $R->run(q` pre = gather(pre, Species, Mixing.Ratio, -Time, -Run) `,
            q` d = rbind(d, pre) `,
    );
}
$R->run(q` d$Species = factor(d$Species, levels = c("O3", "OH", "HO2", "CO", "NO", "NO2")) `);

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, colour = Run)) `,
        q` p = p + geom_line(size = 1) `,
        q` p = p + facet_wrap( ~ Species, scales = "free_y") `,
        q` p = p + theme_tufte() `,
        q` p = p + ggtitle("Mixing Ratio Comparison between Temperatures") `,
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

$R->run(q` CairoPDF(file = "temp_old_vs_new_mixing_ratios.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);
