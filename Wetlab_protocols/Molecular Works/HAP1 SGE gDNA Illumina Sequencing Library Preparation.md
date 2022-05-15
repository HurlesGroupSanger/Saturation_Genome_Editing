# Illumina sequencing library preparation

## Material
gDNA </br>
Primers </br>
KAPA HiFi Ready Mix (Roche, 07958935001) </br>
Qiagen PCR Purification Kit (28104) or 96 well PCR Purification Kit (28181) </br>
ExoI (NEB, M0293L) </br>
Ampure XP beads (Beckman Coulter, A63881) </br>
Qubit dsDNA High Sense kit (Invitrogen, Q32854) </br>
NaOAC (Sodium acetate) (Invitrogen, AM9740) </br>
KAPA Library Quantification Kit (Roche) </br>

## ILL0
PCR from genomic DNA (gDNA)

Primer F: (short primer, 18-30bp, target-specific sequence) </br>
Primer R: (short primer, 18-30bp, target-specific sequence)

### Reaction mix 
100μL of KAPA HiFi Ready Mix </br>
6μL of 10uM Primer F </br>
6μL of 10uM Primer R </br>
4 ug gDNA </br>
Add water to 200μL </br>

### ILL0 PCR conditions
95C: 3min </br>
98C: 20 sec </br>
xxC**: 15 sec </br>
72C: 1.30min </br>
Back to (2) for yy** cycles in total </br>
72C: 1.30min </br>
4C: Forever </br>

** The annealing temperature and cycle number have been pre-optimised by qPCR

### ILL0 and PCR purification protocol
1. Prepare 200μL reactions for each gDNA sample.
2. Split to four 50μL reactions in a 96well PCR plate.
3. Run the PCR with the pre-optimised conditions.
4. Pool the 4 reactions together.
5. Purify PCR products using Qiagen PCR Purification kit (column or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 10μL of NaOAc per 10μL PCR reaction) </br>
  5a. For column purification, use 2 columns per sample. Elude each column in 40μL EB buffer. </br>
  5b. For 96 well plate purification, elude in 80μL EB buffer. </br>
6. Digest the purified DNA with ExoI (NEB, M0293L). 37C for 20 min then 80C for 20 min.  </br>
  6a. For 40μL elutions, make a 30μL reaction (3μL ExoI buffer, 1μL ExoI, 26μL purified DNA) </br>
  6b. For 80μL elutions, make a 60μL reaction (6μL ExoI buffer, 2μL ExoI, 52μL purified DNA) </br>
7. Purify ExoI digested DNA using Qiagen PCR Purification kit (column or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 10μL of NaOAc per 100μL ExoI digested DNA) </br>
  7a. For column purification, use 1 column per sample. Elude each column in 40μL EB buffer. </br>
8. Measure the DNA concentration using Nanodrop. 


## ILL1

Primer F: TCGGCATTCCTGCTGAACCGCTCTTCCGATCT[Target-specific sequence] (IDT, PAGE purified) </br>
Primer R: ACACTCTTTCCCTACACGACGCTCTTCCGATCT[Target-specific sequence] (IDT, PAGE purified)

### Reaction mix
100μL of KAPA HiFi Ready Mix </br>
6μL of 10uM Primer F </br>
6μL of 10uM Primer R </br>
100ng of ILL0 PCR product </br>
Add water to 200μL </br>

### ILL1 PCR conditions
95C: 3min </br>
98C: 20 sec </br>
xxC**: 15 sec </br>
72C: 30 sec </br>
Back to (2) for 6-10 cycles** in total </br>
72C: 30 sec </br>
4C: Forever </br>

** The annealing temperature and cycle number have been pre-optimised by qPCR

### ILL1 and PCR purification protocol
1. Prepare a 200μL reaction for each ILL0 sample.
2. Split to four 50μL reactions in a 96well PCR plate.
3. Run the PCR with the pre-optimised conditions.
4. Pool the 4 reactions together.
5. Purify PCR products using Qiagen PCR Purification kit (column or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 10μL of NaOAc per 100μL PCR reaction) </br>
  5a. For column purification, use 1 column per sample. Elude the column in 40μL EB buffer. </br>
  5b. For 96 well plate purification, elude in 80μL EB buffer. </br>
6. Digest the purified DNA with ExoI (NEB, M0293L). 37C for 20 min then 80C for 20 min.  </br>
  6a. For 40μL elutions, make a 30μL reaction (3μL ExoI buffer, 1μL ExoI, 26μL purified DNA) </br>
  6b. For 80μL elutions, make a 60μL reaction (6μL ExoI buffer, 2μL ExoI, 52μL purified DNA) </br>
7. Purify ExoI digested DNA using Qiagen PCR Purification kit (column or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 10μL of NaOAc per 100μL ExoI digested DNA) </br>
  7a. For column purification, use 1 column per sample. Elude each column in 12μL EB buffer. </br>
8. Measure the DNA concentration using Qubit High Sense dsDNA kit. 

## ILL2

Primer F: 11-base index primer </br>
Primer R: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT (PE 1.0)

### Reaction mix (50uL)
25μL of KAPA HiFi Ready Mix </br>
1.5μL of 10uM Primer F </br>
1.5μL of 10uM Primer R </br>
25ng of ILL1 PCR product </br>
Add water to 50μL

### ILL2 PCR conditions 
95C: 3min </br>
98C: 20 sec </br>
59C: 15 sec </br>
72C: 30sec </br>
Back to (2) for 7 cycles in total </br>
72C: 30sec </br>
4C: Forever

### ILL2 protocol
1. Prepare a 50μL reaction for each sample in a 96 well PCR plate. </br>
2. Run the PCR with the above conditions. </br>
3. Purify the DNA using Ampure beads directly in the 96 well plate using the protocol below.

### Ampure bead purification
1. Warm the Ampure XP bead to room temp. </br>
2. Vortex Ampure XP bead to resuspend. </br>
3. Add 45μL resuspended beads to the 50uL PCR reaction. Mix well by pipetting up and down at least 10 times. Be careful to expel all the liquid out of the tip during the last mix. </br>
4. Incubate samples on bench top for 15 minutes at room temperature. </br>
5. Place the plate on a magnetic stand to separate the beads from the supernatant. </br>
6. After 5 minutes, carefully remove 30μL of the supernatant from the top. The removed supernatant should look clear. Seal the plate and spin at 1000RPM for 20 sec. Most of the bead should be at the bottom and the solution should look clear. Removed the seal carefully and place the plate on the magnetic rack. </br>
7. After 5 minutes, remove all solution. Seal and quick spin the plate. Place the plate on the magnetic rack again.  Remove all remaining solution in the wells.</br>
8. Add 200μl of freshly prepared 80% ethanol to the plate while in the magnetic stand. Incubate at room temperature for 30 seconds, and then carefully remove and discard the supernatant. Be careful not to disturb the beads that contain DNA targets. </br>
9. Repeat Step 8 once for a total of two washes. Seal the plate and spin at 1000RPM for 20 sec. Remove the seal carefully and place the plate on the magnetic rack. Remove traces of ethanol with a p10 pipette tip. Spin the plate again if necessary. </br>
10. Air dry the beads for up to 5 minutes while the plate is on the magnetic stand. Start counting the 5 min after the first plate spinning step from Step 8. </br>
Caution: Do not over-dry the beads. This may result in lower recovery of DNA target. Elute the samples when the beads are still dark brown and glossy looking, but when all visible liquid has evaporated. When the beads turn lighter brown and start to crack they are too dry. </br>
11.	Remove the plate from the magnetic stand. Elute the DNA target from the beads by adding 33μl of Qiagen EB buffer. </br>
12.	Mix well by pipetting up and down 10 times. Incubate for 10 minutes at room temperature. </br>
13.	Spin the plate at 2000RPM for 1 min. </br>
14.	Place the plate on the magnetic stand. After 5 minutes (or when the solution is clear), transfer 30μl to a new plate (ensure that you do not take out any beads). </br>
15.	Measure the DNA concentration using Qubit High Sense dsDNA kit. 

## Pooling and quantification
1. For each library, prepare 15μL of 10ng/uL library. Use Qiagen EB buffer to do the dilution. </br>
2. Use 8μL of the diluted library for 2% gel electrophoresis. 120V, 45min, TAE. You should see a single band at ~400bp. Each library should have similar intensity if the normalisation was done correctly. </br>
3. If the libraries look good on the gel, pool 5μL of each diluted library into one tube. (10ng/uL pooled library). </br>
4. Quantify the pooled library using KAPA Quant qPCR kit.
