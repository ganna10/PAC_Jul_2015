#! /usr/bin/env perl
# Plot tagged and non-tagged species
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use MECCA;
use Statistics::R;

my $base = "/work/users/jco/New_Tagging/MOZART-4_VOC_tagged";
my $mecca = MECCA->new("$base/boxmodel");

my $spc_file = "gas.spc";
my @species = qw(CO HO2 O3 NO2);
my %data;

foreach my $species (@species) {
    my @tagged_species = get_tagged_species($base, $species);
    my $mixing_ratio = $mecca->tracer($species);
    $data{$species}{"Non-tagged"} = $mixing_ratio;
    foreach my $spc (@tagged_species) {
        my $mixing_ratio = $mecca->tracer($spc);
        if ($spc =~ /_/) {
            #$data{$species}{"Ox-tagged"} += $mixing_ratio;
            #} elsif ($spc =~ /_/) {
            $data{$species}{"Tagged"} += $mixing_ratio;
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
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;

$R->set('Time', [  map { $_ } $times->dog ]);
$R->run(q` data = data.frame() `);
foreach my $species (sort keys %data) {
    $R->set('species', $species);
    $R->run(q` pre = data.frame(Time, Species = rep(species, length(Time))) `);
    foreach my $item (sort keys %{$data{$species}}) {
        $R->set('item', $item);
        $R->set('mixing.ratio', [ map { $_ } $data{$species}{$item}->dog ]);
        $R->run(q` pre[item] = mixing.ratio `);
    }
    $R->run(q` pre = gather(pre, Item, Mixing.Ratio, -Time, -Species) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(pre) `);
#print $p, "\n";
$R->set('title', "Real vs Tagged Mixing Ratios");

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Item, group = Item)) `,
        q` plot = plot + geom_line(size = 1) `,
        q` plot = plot + facet_wrap( ~ Species, scales = "free_y") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + ggtitle(title) `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + ylab("Mixing Ratio") `,
        q` plot = plot + xlab("Time (Days)") `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(strip.text = element_text(size = 15, face = "bold")) `,
        q` plot = plot + theme(plot.title = element_text(size = 18, face = "bold")) `,
        q` plot = plot + theme(axis.title = element_text(size = 14, face = "bold")) `,
        q` plot = plot + theme(axis.text = element_text(size = 13)) `,
        q` plot = plot + theme(legend.text = element_text(size = 13)) `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.position = "top") `,
);

$R->set('filename', "tagged_non-tagged_mixing_ratios.pdf");
$R->run(q` CairoPDF(file = filename, width = 10, height = 7) `,
        q` print(plot) `,
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
