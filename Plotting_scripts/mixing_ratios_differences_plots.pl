#! /usr/bin/env perl
# plot time series of argument species mixing ratios in each specification, each mechanism and each model run condition
# Version 0: Jane Coates 24/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $species = $ARGV[0];
die "Need to specify a species\n" unless defined $species;

#my $base = "/work/users/jco/Solvent_Emissions";
my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw(MCM MOZART RADM2);
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
my @runs = qw( Solvents_Only );
my %data;

my $mecca = MECCA->new("$base/RADM2/DE94_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times =  $times(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my $dir = "$base/$mechanism/${speciation}_$run";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $lookup;
            if ($mechanism eq "MOZART" and $species eq "HCHO") {
                $lookup = "CH2O";
            } else {
                $lookup = $species;
            }
            my $mr = $mecca->tracer($lookup); 
            $data{$run}{$mechanism}{$speciation} = $mr(1:$ntime-2) * 1e9;
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(scales) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(grid) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` my.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
foreach my $run (sort keys %data) {
    $R->run(q` data = data.frame() `);
    foreach my $mechanism (sort keys %{$data{$run}}) {
        $R->run(q` pre = data.frame(Time) `);
        foreach my $speciation (sort keys %{$data{$run}{$mechanism}}) {
            $R->set('speciation', $speciation);
            $R->set('mixing.ratio', [ map { $_ } $data{$run}{$mechanism}{$speciation}->dog ]);
            $R->run(q` pre[speciation] = mixing.ratio `);
        }
        $R->run(q` pre = gather(pre, Speciation, Mixing.Ratio, -Time) `);
        $R->set('mechanism', $mechanism) ;
        $R->run(q` pre$Mechanism = rep(mechanism, length(pre$Mixing.Ratio)) `,
                q` data = rbind(data, pre) `,
        );
    }
    $R->set('filename', "${run}_${species}_mixing_ratios.pdf");
    $R->set('title', "$species Mixing Ratios");
    $R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
    $R->run(q` data = mutate(data, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2"))) `);
    
    $R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Speciation, group = Speciation)) `,
            q` plot = plot + geom_line(size = 1) `,
            q` plot = plot + facet_wrap( ~ mechanism) `,
            q` plot = plot + theme_tufte() `,
            q` plot = plot + ggtitle(title) `,
            #q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
            q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
            q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
            q` plot = plot + xlab("Time (Days)") `,
            q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
            q` plot = plot + theme(strip.text = element_text(size = 15, face = "bold")) `,
            q` plot = plot + theme(plot.title = element_text(size = 18, face = "bold")) `,
            q` plot = plot + theme(axis.title = element_text(size = 14, face = "bold")) `,
            q` plot = plot + theme(axis.text = element_text(size = 13)) `,
            q` plot = plot + theme(legend.text = element_text(size = 13)) `,
            q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
            q` plot = plot + theme(legend.title = element_blank()) `,
            #q` plot = plot + theme(legend.position = c(1.03, 1.03)) `,
            #q` plot = plot + theme(legend.justification = c(1, 1)) `,
            q` plot = plot + scale_colour_manual(values = my.colours) `,
            #q` plot = plot + guide(guides = guide_legend(position = "horizontal")) `,
            q` plot = plot + theme(legend.position = "top") `,
    );

    $R->run(q` CairoPDF(file = filename, width = 12, height = 8) `,
            q` print(plot) `,
            q` dev.off() `,
    );
};

$R->stop();
