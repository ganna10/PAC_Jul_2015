#! /usr/bin/env perl
# Time series of all tagged O3 mixing ratios stacked to give total O3
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/New_Tagging";
my @mechanisms = qw( MOZART-4 );
my $species = "O3";
my %data;

my %category_mapping = (
    MOZART  =>  {   BIGALK  => [ '0.285 NC4H10', '0.151 IC4H10', '0.146 NC5H12', '0.340 IC5H12', '0.048 NC6H14', '0.020 NC7H16', '0.010 NC8H18' ],
                    BIGENE  => [ '0.333 BUT1ENE', '0.667 MEPROPENE' ],
                    TOLUENE => [ '0.166 BENZENE', '0.478 TOLUENE_M', '0.073 EBENZ', '0.142 MXYL', '0.069 OXYL', '0.073 PXYL' ],
                }
);

my $mecca = MECCA->new("$base/MOZART-4_VOC_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;
$time = $time(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    my $dir = "$base/${mechanism}_VOC_tagged";
    my $mecca = MECCA->new("$dir/boxmodel");
    my @tagged_species = get_tagged_species($dir, $species);
    foreach my $spc (@tagged_species) {
        my ($base, $tag) = split /_X_/, $spc;
        my $mixing_ratio = $mecca->tracer($spc);
        $mixing_ratio = $mixing_ratio(1:$ntime-2);
        $data{$mechanism}{$tag} = $mixing_ratio;
    }
}

foreach my $mechanism (keys %data) { #allocate MOZART species to Functional Groups
    foreach my $spc (keys %{$data{$mechanism}}) {
        if ($spc eq "C2H6" or $spc eq "C3H8" or $spc eq "BIGALK") {
            $data{$mechanism}{"Alkanes"} += $data{$mechanism}{$spc};
            delete $data{$mechanism}{$spc};
        } elsif ($spc eq "C2H4" or $spc eq "C3H6" or $spc eq "BIGENE") {
            $data{$mechanism}{"Alkenes"} += $data{$mechanism}{$spc};
            delete $data{$mechanism}{$spc};
        } elsif ($spc eq "TOLUENE") {
            $data{$mechanism}{"Aromatics"} += $data{$mechanism}{$spc};
            delete $data{$mechanism}{$spc};
        } elsif ($spc eq "ISOP") {
            $data{$mechanism}{"Isoprene"} += $data{$mechanism}{$spc};
            delete $data{$mechanism}{$spc};
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
);

$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` d = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time) `);
    foreach my $voc (sort keys %{$data{$mechanism}}) {
        $R->set('voc', $voc);
        $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$voc}->dog ]);
        $R->run(q` pre[voc] = mixing.ratio  * 1e9 `); #convert to ppb
    }
    $R->run(q` pre = gather(pre, VOC, Mixing.Ratio, -Time) `,
            q` d = rbind(d, pre) `,
    );
}

#$R->run(q` d$VOC = factor(d$VOC, levels = c("INI", "XTR", "CO", "Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Propene", "Butene", "2-Methylpropene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Ethylbenzene")) `,
#        q` my.colours = c( "Ethane" = "#696537", "Propane" = "#f9c500", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb", "Ethene" = "#86b650", "Propene" = "#6c254f", "Butene" = "#ee6738", "2-Methylpropene" = "#58691b", "Isoprene" = "#8ed6d5", "Benzene" = "#1c3e3d", "Toluene" = "#c65d6c", "m-Xylene" = "#888a87", "o-Xylene" = "#0e5c28", "p-Xylene" = "#b569b3", "Ethylbenzene" = "#2c9def", "Methane" = "#000000", "INI" = "#c9a415", "XTR" = "#f3aa7f", "CO" = "#d94a80" ) `,
#);
$R->run(q` d$VOC = factor(d$VOC, levels = c("INI", "XTR", "CO", "CH4", "Alkanes", "Alkenes", "Isoprene", "Aromatics")) `,
        q` my.colours = c("INI" = "#6c254f", "XTR" = "#f9c500", "CO" = "#0e5c28", "CH4" = "#2b9eb3", "Alkanes" = "#ef6638", "Alkenes" = "#86c650", "Isoprene" = "#8c1531", "Aromatics" = "#0352cb") `,
);

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, fill = VOC, order = VOC)) `,
        q` p = p + geom_area(position = "stack") `,
        q` p = p + geom_area(position = "stack", colour = "black", show_guide = FALSE) `,
        q` p = p + theme_tufte() `,
        q` p = p + ggtitle("Source Allocation of Ozone Mixing Ratios") `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` p = p + scale_y_continuous(expand = c(0, 1e-9)) `,
        q` p = p + ylab("Mixing Ratio (ppbv)") `,
        q` p = p + xlab("Time (Days)") `,
        q` p = p + scale_fill_manual(values = my.colours, limits = rev(levels(d$VOC))) `, 
        q` p = p + theme(plot.title = element_text(size = 18, face = "bold")) `,
        q` p = p + theme(axis.title = element_text(size = 14, face = "bold")) `,
        q` p = p + theme(axis.text = element_text(size = 13)) `,
        q` p = p + theme(legend.text = element_text(size = 13)) `,
        q` p = p + theme(legend.position = "top") `,
);

$R->run(q` CairoPDF(file = "O3_mixing_ratio_components.pdf", width = 10, height = 7) `,
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

sub get_contributions {
    my ($spc, $mechanism, $mixing_ratio) = @_;
    $mechanism = "MOZART" if ($mechanism =~ /MOZ/);

    my %data;
    if (exists $category_mapping{$mechanism}) { 
        if (exists $category_mapping{$mechanism}{$spc}) {
            my $list = $category_mapping{$mechanism}{$spc}; 
            foreach (@$list) {
                my ($fraction, $VOC) = split /\s/, $_;
                $VOC = get_voc_name($VOC);
                $data{$VOC} = $fraction * $mixing_ratio;
            }
        } else {
            my $VOC = get_voc_name($spc);
            $data{$VOC} = $mixing_ratio;
        }
    }
    return \%data;
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
