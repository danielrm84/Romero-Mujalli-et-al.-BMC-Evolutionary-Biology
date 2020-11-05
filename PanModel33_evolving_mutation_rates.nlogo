;**************************************************************************************
;          Project Name: Evolutionary adaptive responses to rapid environmental change
;                                  UNIVERSITY OF POTSDAM
;
; Author:          Daniel Romero Mujalli
;
; Written:         27. 6.2016
; Last update:     31. 8.2018
;
; Type of model:   Individual-based model (IBM)
;
; Summary:         The purpose of the model is to investigate whether different
;              mutation rates and effect size of mutations may arise depending on
;              the scenario of environmental change, and on the type of organism.
;
;              REMAINING TASKS:
;             - mutation rate per locus involve the two alleles? e.g., if mu = 1
;               both alleles mutate? if yes, need to modify evolution of m-rate!
;
;              NOTES / COMMENTS / QUESTIONS:
;             - values of strength of selection of 2.5; 10 and 20 are the ones that
;               yield values similar to those reported by Bjoerklund2009 when the
;               phenotype departs from optimum in one phenotypic standard deviation
;
;**************************************************************************************
;**************************************************************************************
;                 Copyright © 2017 Daniel Romero Mujalli
;
;        This program is free software: you can redistribute it and/or modify
;        it under the terms of the GNU General Public License as published by
;        the Free Software Foundation, either version 3 of the License, or
;        (at your option) any later version.
;
;        This program is distributed in the hope that it will be useful,
;        but WITHOUT ANY WARRANTY; without even the implied warranty of
;        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;        GNU General Public License for more details.
;
;        You should have received a copy of the GNU General Public License
;        along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;**************************************************************************************

; load r-extension
extensions [r]

;global parameters of the model:
globals
[
  genetic-mean          ; initial genetic mean of the population
  env-effect-mean       ; mean random environmental effect
  initial-env-optimum   ; the initial phenotypic optimum as given by the environment
  size-of-environment   ; the preferred size of the environment (distance from the center)
  female-male-ratio     ; controls the sex ratio in the population. Example: a value of 0.7
                        ; means 70% females and 30% males in the population
  env-effect-variance   ; variance of environmental effect
  extinction?           ; boolean: true if the population is extinct
  pop-gv                ; population level additive genetic variance
  density-compensation  ; governs the strength of density dependence effect
  strength-selection    ; strength of selection used in fitness function
  ;time-limit            ; desired time limit (in generations)
]

; properties of the individuals
turtles-own
[
  phenotype            ; phenotype 'z'
  genetic-component    ; genetic component 'a' of the phenotype
  environmental-effect ; random environmental effect 'e'
  fitness              ; fitness wi of individual i
  fecundity            ; fecundity lambda of individual i
  stage                ; whether "adult" or "juvenile", as string
  sex                  ; sex: "female" or "male"
  reproduced?          ; boolean.
  dna-strain1          ; first strain of the chromosome
  dna-strain2          ; second strain of the chromosome
  my-mut-rate          ; mutation rate of this turtle
  my-mut-effect-size   ; mutation effect size of this turtle
  mu-dist-effect-size  ; mean of the distribution of mutation effects
  MR                   ; just to report the MR
]

; properties of the environment (patches)
patches-own
[
  optimum              ; phenotypic optimum tita as given by the environment
  noise                ; variance of the environmental optimum
  degree-maladaptation ; degree of maladaptation
  time-to-extinction   ; store the time when the population went extinct
]






;**************************************************************************************
;        TO SETUP
;**************************************************************************************
; This procedure set the values for the global parameters of the model, the optimum
; phenotype as given by the environment (here the environment consists of a group of
; patches with green color. Then each individual is created with own values for the
; genetic and environmental components. These two components are used to calculate
; the phenotye. Finally, the distribution of phenotypes of the population is plotted.
to setup

  ; reset setup
  clear-all
  r:clear
  reset-ticks

  ; set values for global parameters
  set-global-parameters

  ; set the environment:
  ; patch at the center (0,0) and those in radius 10 from the center are use as
  ; the environment to test for local adaptation
  ; more than one patch were selected for better visualization of the world
  ask patch 0 0
  [
    set pcolor green
    set optimum initial-env-optimum
    set noise 0
    set degree-maladaptation compute-degree-maladaptation
    set time-to-extinction 0
    ask patches with [distance myself <= size-of-environment]
    [
      set pcolor green
      set optimum initial-env-optimum
      set noise 0
      set degree-maladaptation compute-degree-maladaptation
      set time-to-extinction 0
    ]
  ]

  ; create a population of N individuals,
  ; initialize individuals (the turtles) and
  ; calculate phenotype (update phenotype)
  create-turtles population-size
  [
    set color blue
    move-to one-of patches with [pcolor = green]
  ]

  ; Set initial mutation rate and mutation effect size
  ask turtles [ set-initial-mut-rate-and-effect-size]

  ask turtles [ set MR 10 ^ (-1 * mean my-mut-rate) ]

  ; set initial conditions according to the method for simularing genetics
  if-else (how-genetics? = "implicit")
  [ ask turtles [ set-initial-conditions-implicit-genetics ] ]
  [
    if-else (how-genetics? = "explicit")
    [ ask turtles [ set-initial-conditions-explicit-genetics ] ]
    [ error "undefined method for modelling genetics" stop ] ;catch exception
  ]
  ; update phenotype
  ask turtles [ update-phenotype ]

  ; update mean of the distribution of effect size
  ask turtles [ update-mean-dist-effects ]

  ; calculate fitness of each individual
  if-else( fitness-function = "Bjoerklund2009")
  [
    ask turtles [ check-fitness-Bjoerklund2009 ]
    ; DEBUG:
    ;print "fitness according to Bjoerklund"
  ]
  [ ; else, negative-exponential
    if-else (fitness-function = "negative-exponential")
    [ ask turtles [ check-fitness-negative-exponential ] ]
    [ error "fitness function not identified" stop ] ; else, error message
  ]

  ; update fitness-plot
  plot-mean-fitness

  ;update-output:
  update-output

end






;**************************************************************************************
;        TO GO
;**************************************************************************************
to go

  ; iterations are counted at the beginning of the loop.
  ; Each iteration means one generation. The model
  ; uses non-overlapping generations:
  tick

  ; Environmental scenario:
  ; update-environment: updates optimum phenotype as given by the-environment according
  ; to the scenario of environmental change.
  ; If specified, the environmental change starts after a given number of
  ; iterations/generations (mutation selection balance/mutation-selection drift balance)
  ;if ( ticks > time-to-balance )
  ;[
    ;ask patches with [ pcolor = green ] [ update-environment ]
    ask patch 0 0 [ update-environment ]
  ;]

  ; update-stage, reproductive status and sex of individuals:
  ask turtles [   update-stage-sex  ]

  ; calculate fitness of each individual
  if-else( fitness-function = "Bjoerklund2009")
  [
    ask turtles [ check-fitness-Bjoerklund2009 ]
    ; DEBUG:
    ;print "fitness according to Bjoerklund"
  ]
  [ ; else, negative-exponential
    if-else (fitness-function = "negative-exponential")
    [ ask turtles [ check-fitness-negative-exponential ] ]
    [ error "fitness function not identified" stop ] ; else, error message
  ]

  ; update fitness-plot
  plot-mean-fitness

  ; calcualte fecundity of each individual
  ask turtles [ set-fecundity-Bjoerklund ]

  ; reproduction
  ; store population-level additive genetic variance
    if (how-genetic-variance = "population-level")
    [
      if-else ( count turtles > 1)
      [
        set pop-gv variance [genetic-component] of turtles with [stage = "adult"]
      ]
      [ set pop-gv 0 ] ; else, if only one turtle
    ]
  ; lottery polygyny: females reproduce only once, but males can repeat
  ask turtles with [stage = "adult" and sex = "female" ]
  [
    if-else (how-genetics? = "implicit")  [ reproduce-implicit-genetics ]
    [ if-else (how-genetics? = "explicit")[ reproduce-explicit-genetics ]
      [ error "undefined method of how-genetics?" stop ] ; catch exception
    ]
  ]

  ask turtles [ set MR 10 ^ (-1 * mean my-mut-rate) ]

  ; mortality of adults
  ask turtles with [ stage = "adult" ] [ die ]

  ; update-phenotype
  ask turtles [ update-phenotype ]

  ; update mean of the distribution of effect size
  if ( beneficial-mutations != 0.5 )
  [
    if (evolving-mut-rate? = "true" or
        evolving-mut-rate? = "only-mut-effect-size")
    [ ask turtles [ update-mean-dist-effects ] ]
  ]

  ; compute the degree of maladaptation according to Björklund et al (2009)
  ask patch 0 0
  [
    ; DEBUG:
    ;write "degree-maladaptation before: " print degree-maladaptation
    let new-degree-maladaptation compute-degree-maladaptation

    set degree-maladaptation new-degree-maladaptation

    ; DEBUG:
    ;write "degree-maladaptation after: " print degree-maladaptation

    ask patches with [distance myself <= size-of-environment]
    [ set degree-maladaptation new-degree-maladaptation ]
  ]

  ; check time-to-extinction:
  if (count turtles < 1 and extinction? = false)
  [
    ask patches with [ pcolor = green ] [ set time-to-extinction (ticks - time-to-balance) ]
    set extinction? true
  ]

  ; update-output:
  update-output

  ; end simulation:
  if (extinction? = true or (ticks - time-to-balance) > time-limit) [ stop ]

end






;**************************************************************************************
;        TO SET-GLOBAL-PARAMETERS
;**************************************************************************************
to set-global-parameters

  ;set time-limit           30 ; in generations
  set genetic-mean         0.0
  set env-effect-mean      0.0
  set initial-env-optimum  0.0
  set size-of-environment 10
  set female-male-ratio    0.5 ; 0.5 means random set of sex
  set extinction? false
  ; set the variance of the random environmental effect "ve" on the development of the
  ; trait. If how-heritability? = fixed-value, then ve is computed according to given h2 and gv
  ; as in Bjoerklund et al (2009).
  ; Otherwise, if h2 emerge, then ve is assumed constant, ve = 1 as in Vincenzi (2014)
  if-else (how-heritability? = "fixed-value")
  [ set env-effect-variance   compute-env-effect-variance (heritability) (genetic-variance)]
  [ set env-effect-variance 1 ]

  ; set level of density compensation, which is simulated as in Bjoerklund et al (2009),
  ; according to the density dependence effeect selected by the user: weak, normal, strong,
  ; or very strong
  ; Density compensation allows to simulate different population dynamics. For example,
  ; stable or oscillatory population dynamics.
  if-else (density-dependence-effect = "weak")
  [ set density-compensation 0.5 ]
  [;else
    if-else (density-dependence-effect = "normal")
    [ set density-compensation 1 ]
    [;else
      if-else (density-dependence-effect = "strong")
      [ set density-compensation 1.8]
      [;else
        if-else (density-dependence-effect = "very strong")
        [ set density-compensation 2.5 ]
        [ if (density-dependence-effect = 0)
          [ error "Warning! no density dependence!" stop ]
        ]
      ]
    ]
  ]
  ; DEBUG:
  ;write "density dependence effect: " print density-dependence-effect
  ;write "value of density compensation: " print density-compensation


  ; set the strength of selection according to fitness function and type of organism
  ; whether specialist, moderate, generalist
  if-else (fitness-function = "Bjoerklund2009")
  [
    ; values are set as in Bjoerklund et al (2009)
    if-else (type-organism = "specialist") [ set strength-selection 10 ] ;10 ] 2.5
    [ ;else
      if-else (type-organism = "moderate") [ set strength-selection 20 ] ;20 ] 10
      [;else
        if-else (type-organism = "generalist") [ set strength-selection 40 ] ;40 ] 20
        [;else, catch exception
          error "Unidentified type of organism ..." stop
        ]
      ]
    ]
  ]
  [;else: exponential function
    ; values are set as in Burger and Lynch (1995): strength of selection
    ; gamma ^ 2 = 1; 9; 99. Instead of 99, I am using gamma = 9 ^ 2 = 81
    if-else (fitness-function = "negative-exponential")
    [
      if-else (type-organism = "specialist") [ set strength-selection 1 ] ; 1
      [ ;else
        if-else (type-organism = "moderate") [ set strength-selection 2.2 ] ; 3
        [;else
          if-else (type-organism = "generalist") [ set strength-selection 3.2 ] ; 9
          [;else, catch exception
            error "Unidentified type of organism ..." stop
          ]
        ]
      ]
    ]
    [;else; catch exception
      error "Unidentified fitness function" stop
    ]
  ]

  ; DEBUG
  ;write "strength of selection: " print strength-selection
  ;write "env-effect-variance: " print env-effect-variance

  ; end of parameter values

end






;**************************************************************************************
;        TO SET-INITIAL-MUT-RATE-AND-EFFECT-SIZE
;**************************************************************************************
to set-initial-mut-rate-and-effect-size

  ; mutation rate:
  if-else (how-mut-rate? = "input-value" or evolving-mut-rate? = "false")
  [ set my-mut-rate list (mut-rate-per-locus) (mut-rate-per-locus) ]
  [ ;else, random initial value
    if-else (method-for-MR = "decimal")
    [
      set my-mut-rate n-values 2 [((random 9) + 1) * 10 ^ (-1 * ((random 4) + 1)) ]
    ]
    [
      if-else (method-for-MR = "exponent")
      ;[ set my-mut-rate (list 3 3)]
      [ set my-mut-rate n-values 2 [ (random 3) + 2 ] ]
      ;[ set my-mut-rate n-values 2 [ (random-float 3) + 2 ] ] ; as in Cobben et al 2017
      [ error "unspecified method for mutation rate MR" stop]
    ]
  ]

 ; mutation effect size:
  if-else (how-mut-effect-size? = "input-value" or evolving-mut-rate? = "false")
  [ set my-mut-effect-size list (mut-effect-size) (mut-effect-size) ]
  [ ;else, random initial value
    set my-mut-effect-size n-values 2 [((random 9) + 1) * 10 ^ (-1 * (( random 4) + 1)) ]
  ]

  ;DEBUG:
  ;write "mut-rate: " print mean my-mut-rate
  ;write "effect-size: " print mean my-mut-effect-size

end







;**************************************************************************************
;        TO SET-INITIAL-CONDITIONS-IMPLICIT-GENETICS
;**************************************************************************************
; This function set the initial conditions for both, the genetic and environmental
; components of the phenotype for each individual (turtle).
; Here the genetics is implicitly modelled.
; Important!: random-normal function uses the standard deviation = sqrt(variance) of the
; distribution.
; The function works in a turtle context, example: ask turtles [ set-initial-cond... ]
to set-initial-conditions-implicit-genetics

  ; genetic component 'a':
  set genetic-component random-normal genetic-mean sqrt (genetic-variance)

  ; environmental component 'e':
  set environmental-effect random-normal env-effect-mean sqrt (env-effect-variance)

  ; DEBUG:
  ;write "turtle_" write who print ": "
  ;write "genetic component: " print genetic-component
  ;write "environ component: " print environmental-effect

end







;**************************************************************************************
;        TO SET-INITIAL-CONDITIONS-EXPLICIT-GENETICS
;**************************************************************************************
; this function sets the initial conditions for the genetic and environmental component
; of each individual turtle.
; The genetics is modelled explicitly, considering the number of loci.
; Three methods for the initialization of allele values were implemented. Only one must
; be active.
; alleles can be either:
; - real allele values set according to normal distribution as in Vincenzi (2014)
; - real allele values randomly set
; - random set of ones and zeros (on / off)
; (loci effect on phenotype is assumed to be additive).
; This function works in a turtle context, example: ask turtles [set-initial-condi... ]
to set-initial-conditions-explicit-genetics

  ; set allele values for each locus. Notice that the two dna strains are simulated
  ; separately

  ; population is assumed locally adapted. Allele values from a normal distribution
  ; as in Vincenzi (2014)
  ; N(0, va), where mean = 0 = initial environmental optimum, and
  ; va is the additive genetic variance per locus at the start of the simulation, given by
  ; VG (genetic-variance) / 2*L, L is the number-of-loci, and VG is the additive genetic
  ; variance of the quantitative trait at the start of the simulation
  ; (THIS METHOD is based on the continuum-of-alleles model of Kimura 1970: an introduction
  ; to population genetics theory (cited in Vincenzi 2014) )
  ; Vincenzi (2014) simulate additive effect through the addition of allele values
  set dna-strain1 n-values number-of-loci [ random-normal 0
                                           sqrt (genetic-variance / (2 * number-of-loci)) ]
  set dna-strain2 n-values number-of-loci [ random-normal 0
                                           sqrt (genetic-variance / (2 * number-of-loci)) ]

  ; alleles take values from uniform distribution in range (-1. 1)
  ;set dna-strain1 n-values number-of-loci [ (random-float 2) - 1 ]
  ;set dna-strain2 n-values number-of-loci [ (random-float 2) - 1 ]

  ; another method: alleles take values of 1 or 0 (on / off)
  ; this method might not be appropriated for small number of loci
  ; set dna-strain1 n-values number-of-loci [ random 2 ]
  ; set dna-strain2 n-values number-of-loci [ random 2 ]

  ; DEBUG
  ;type "dna-strain1: " print dna-strain1
  ;type "dna-strain2: " print dna-strain2
  ;show item (number-of-loci - 1) dna-strain1

  ; set the value for the genetic component 'a' according to the genome
  ; loci effect is assumed to be additive
  set genetic-component (sum dna-strain1) + (sum dna-strain2)

  ; set environmental component 'e':
  set environmental-effect random-normal env-effect-mean sqrt (env-effect-variance)

  ; DEBUG:
  ;write "turtle_" write who print ": "
  ;write "genetic component: " print genetic-component
  ;write "environ component: " print environmental-effect
end






;**************************************************************************************
;        TO UPDATE-PHENOTYPE
;**************************************************************************************
; This function calculates (updates) the phenotype 'z' of the turtle following
; z = a + e
; where 'a' is the genetic-component and 'e' the environmental-effect in development
; The function works in a turtle context, example: ask turtles [ update-phenotype ]
to update-phenotype

  set phenotype genetic-component + environmental-effect

  ;DEBUG:
    ;write "phenotype z of turtle " write who print ":"
    ;print phenotype
    ;write "genetics: " print genetic-component
    ;write "env-effect: " print environmental-effect

end






;**************************************************************************************
;        TO UPDATE-STAGE-SEX
;**************************************************************************************
; update the stage of the individuals to adulthood, and set sex randomly.
; The function works in a turtle context, example: ask turtles [ update-stage-sex ]
to update-stage-sex

  ; set stage to adulthood:
  set stage "adult"
  ; set reproduced? to false
  set reproduced? false
  ; set sex randomly:
  if-else (random-float 1 <= female-male-ratio)
  [ set sex "female" ]
  [ set sex  "male"  ] ; else, male

end






;**************************************************************************************
;        TO CHECK-FITNESS-BJOERKLUND2009
;**************************************************************************************
; Fitness wi of individual i is calculated according to Björklund et al. (2009):
; wi = 1 - [(zi - tita(t))^2]/gamma
; where tita(t) is the optimum phenotype as given by the environment at time t, and
; gamma is the strength of selection
; the maximum fitness is 1
; The function works in a turtle context
to check-fitness-Bjoerklund2009

  if-else([pcolor] of patch-here = green)
  [
    set fitness (1 - (((phenotype - [optimum] of patch-here) ^ 2) / strength-selection))
  ]
  [ error "turtle outside green patches" stop ]

  ; DEBUG
  ;write "fitness of turtle " write who write ": " print fitness
  ;write "optimum of patch-here: " print [optimum] of patch-here
  if (fitness > 1)
  [
    write "phenotype of turtle : " write who write ": " print phenotype
    write "optimum: " print [optimum] of patch-here
    write "gamma: " print strength-selection
    write "fitness: " print fitness
    error "warning! fitness value is greater than 1" stop
  ]

end






;**************************************************************************************
;        TO CHECK-FITNESS-NEGATIVE-EXPONENTIAL
;**************************************************************************************
; Fitness wi of individual i is calculated according to a negative exponential function
; as in Burger & Lynch (1995):
; wi = exp [ - (zi - tita(t))^2]/2*gamma^2
; where tita(t) is the optimum phenotype as given by the environment at time t, and
; gamma is the strength of selection.
; wi is in range (0, 1) => maximum fitness of 1
; The function works in a turtle context
to check-fitness-negative-exponential

  let tita 0

  if-else([pcolor] of patch-here = green)
  [ set tita [optimum] of patch-here ]
  [ error "turtle outside green patches!" stop]

  let gamma strength-selection

  set fitness exp ( - ( (phenotype - tita) ^ 2) / (2 * (gamma) ^ 2) )

  ; DEBUG
  ;write "fitness of turtle " write who write ": " print fitness
  ;write "optimum of patch-here: " print [optimum] of patch-here

  ; catch error:
  if (fitness > 1)
  [
    write "phenotype of turtle : " write who write ": " print phenotype
    write "optimum: " print [optimum] of patch-here
    write "gamma: " print strength-selection
    write "fitness: " print fitness
    error "warning! fitness value is greater than 1" stop
  ]

end






;**************************************************************************************
;        TO SET-FECUNDITY-BJOERKLUND
;**************************************************************************************
; Fecundity lambda is calculated according to Björklund et al. (2009) for each
; individual:
; lambda = wi*exp[alpha(1 - N/K)]
; where wi is the fitness of individual i, alpha is the density-compensation, N the
; population size, K the carrying capacity and exp the exponential function
; The function works in a turtle context
to set-fecundity-bjoerklund

  let alpha density-compensation
  let N count turtles with [ stage = "adult" ]
  let K carrying-capacity
  let wi fitness

  set fecundity wi * ( exp ( alpha * (1 - N / K ) ) )

  ; DEBUG
  ;write "fecundity of turtle " write who write ": " print fecundity
  ;write "population size: " print N

end







;**************************************************************************************
;        TO REPRODUCE-IMPLICIT-GENETICS
;**************************************************************************************
; This function implements sexual reproduction according to Björklund et al. (2009).
; Female individuals randomly select a partner of opposite sex to mate. The fitness of
; a pair is the sum of the fitness value of the two parents. Then fecundity is calculated
; considering density dependence. The number of offspring is drawn from a poisson
; distribution centered on the value of fecundity.
;
; Genetics is implicit according to the infinitesimal model of quantitative genetics which
; assumes that traits are affected by a large number of loci of additive effects.
; Therefore trait inheritance can be approximated using a normal distribution with mean
; parents trait value, and variance: half the genetic variance of the parents. Mutations
; and recombination are implicit.
;
; The function works in a turtle context, example: ask turtles [ reproduce ]
to reproduce-implicit-genetics

  ; store the sex of this individual
  ;let my-sex sex

  ; SEXUAL REPRODUCTION:
  ; pick a random partner of opposite sex:
  let my-partner one-of turtles with [ stage = "adult" and
                                       sex = "male"   ];sex != my-sex
                                       ;reproduced? = false   ]

  ; test whether there is a partner available to reproduce:
  if ( my-partner != nobody )
  [
    ; DEBUG:
    ;write "me "print who
    ;write "partner " print [who] of my-partner
    ;type "my-sex: " print sex
    ;type "partner sex: " print [sex] of my-partner
    ;if(reproduced? = true) [ print "me"]
    ;if([reproduced?] of my-partner = true) [print "partner"]

    ;create a list containing the genetic-component of both parents:
    let list-genetic-parents list (genetic-component) ([genetic-component] of my-partner)

    ;DEBUG:
    ;write "my-genetic-component: " print genetic-component
    ;write "my-partner-genetic-comp: " print [genetic-component] of my-partner
    ;write "list-genetic-components-of-parents: " print list-genetic-parents

    ; calculate the mean parental genetic-component which will be inherited by offspring
    let genetic-mean-of-parents mean ( list-genetic-parents )

    ;DEBUG:
    ;write "my genetic-component: " print genetic-component
    ;write "partner genetic-component: " print [genetic-component] of my-partner
    ;write "genetic-mean: " print genetic-mean-of-parents

    ; calculate genetic variance among parents:
    let genetic-variance-of-parents variance ( list-genetic-parents )
    ; DEBUG:
    ;write "genetic-variance-of-parents: " print genetic-variance-of-parents
    ;type "adults "print count turtles with [stage = "adult"]
    ;type "total " print count turtles
    ; set the genetic-variance for the offspring according to the selected method
    let my-genetic-variance 0

    if-else ( how-genetic-variance = "parental-level")
    [
      set my-genetic-variance (genetic-variance-of-parents);parental-level
                                                           ;(Bjoerklund2009)
    ]
    [; else,
      if-else (how-genetic-variance = "parameter")
      [ set my-genetic-variance genetic-variance ] ; as in Reed et al. (2011):
      ;the population-level additive genetic variance is an input parameter (set by user)
      [; else
        if-else (how-genetic-variance = "population-level")
        [
        ; method as in Vincenzi & Piotti (2014): genetic variance is the total
        ; additive genetic variance for the trait at the population level
          set my-genetic-variance pop-gv
        ]
        ; else, catch exception
        [ error "undefined method for genetic variance" stop ]
      ]
    ]
    ; the genetic variance is simulated as in Vincenzi et al (2012). The infinitesimal
    ; model is modified to account for the decline in additive genetic variance and the
    ; new input of variation through mutation (mutational-variance).
    ; This method also accounts for offspring additive genetic variance equals half the
    ; additive genetic variance
    set my-genetic-variance (1 / 2) * (my-genetic-variance + mutational-variance)

    ; calculate the fecundity of the breeding pair pair-fecundity as in Björklund et al.
    ; Björklund(2009) assumed that the fitness of the pair was the sum of the fitness
    ; values of the two parents (wsum) which is equivalent to the sum of their fecundities.
    ; Given that wsum is used, the population is saturated above the carrying capacity.
    let fecundity-of-pair ( fecundity + [ fecundity ] of my-partner )

    ;DEBUG:
    ;write "my-fecundity: " print fecundity
    ;write "fecundity-of-partner: " print [fecundity] of my-partner
    ;write "fecundity-of-pair: " print fecundity-of-pair

    let number-of-offspring random-poisson fecundity-of-pair
    ; DEBUG
    ;write "nr. offspring: " print number-of-offspring

    if (number-of-offspring >= 1) ; reproduction occurs
    [
      hatch number-of-offspring
      [ ; according to netlogo dictionary, each new turtle inherits of all its
        ; variables, including its location, from its parent, except [who]

        ; set stage as juvenile
        set stage "juvenile"
        ; set random orientation of the newborn and move to a green patch
        set heading random 360
        move-to one-of patches with [pcolor = green]

        set genetic-component random-normal genetic-mean-of-parents ; mean
                                sqrt my-genetic-variance

        ; DEBUG
        ;type "genetic-component of offspring: " print genetic-component
        ;type "genetic-variance: "  print my-genetic-variance

         ; environmental component 'e':
        ; Set the variance of environmental effect ve according to method of heritability
        ; if h2 = fixed, compute the variance of environmental effect according to the
        ; genetic variance and level of heritability
        ; else, ve is a constant parameter
        let ve 0 ; initialization
        if-else (how-heritability? = "fixed-value")
        [ set ve compute-env-effect-variance (heritability) (my-genetic-variance)]
        [ set ve env-effect-variance                                             ]

        set environmental-effect random-normal env-effect-mean sqrt ve

        ; DEBUG:
        ;write "genetic-variance: " print my-genetic-variance
        ;write "env-effect-variance: " print ve
      ] ; end of inheritance

   ] ; end of hatching n offspring ( reproduction)

    ; update reproduced? status of the partner to true
    ask my-partner [ set reproduced? true ]
    ;type "genetic-component of partner " print [genetic-component] of my-partner
  ]; end of if statement: whether a partner is available for reproduction

  ; update reproduced? of this turtle or individual to true:
    set reproduced? true
    ;type "my-genetic-component " print genetic-component

end






;**************************************************************************************
;        TO REPRODUCE-EXPLICIT-GENETICS
;**************************************************************************************
; This function implements sexual reproduction and lottery polygyny.
; Females randomly select a partner of opposite sex to mate. The fitness of a pair
; is the sum of the fitness value of the two parents. Then fecundity is calculated
; considering density dependence. The number of offspring is drawn from a poisson
; distribution centered on the value of fecundity.
;
; Genetics is explicit. This means that loci are explicitly modelled. Loci effect on
; trait value is assumed to be additive.
; Recombination is simulated through the random selection of parental allele values.
; Then, mutations take place with probability mut-rate.
;
; The function works in a turtle context, example: ask turtles [ reproduce ]
to reproduce-explicit-genetics

  ; store the sex of this individual
  ;let my-sex sex
  ; store the loci (dna strains) of this turtle (parent 1)
  let my-loci1 dna-strain1
  let my-loci2 dna-strain2

  ; SEXUAL REPRODUCTION:
  ; pick a random partner of opposite sex:
  let my-partner one-of turtles with [ stage = "adult" and
                                       sex = "male" ]

  ; test whether there is a partner available to reproduce:
  if ( my-partner != nobody )
  [
    ; DEBUG:
    ;write "me "print who
    ;write "partner " print [who] of my-partner
    ;type "my-sex: " print sex
    ;type "partner sex: " print [sex] of my-partner
    ;if(reproduced? = true) [ print "me"]
    ;if([reproduced?] of my-partner = true) [print "partner"]

    ; store loci of my-partner (parent 2)
    let partner-loci1 [dna-strain1] of my-partner
    let partner-loci2 [dna-strain2] of my-partner

    ;DEBUG:
    ;write "Loci parent 1: " write my-loci1 write " and " print my-loci2
    ;write "Loci parent 2: " write partner-loci1 write " and " print partner-loci2

    ; set mut-rate and mut-effect-size for each parent, for the process of inheritance
    ; mutations take place during the creation of gametes. Therefore, the action of the
    ; mutator locus of each parent needs to be considered.
    ; parent 1:
    let mut-rate-parent1 my-mut-rate
    let MV-parent1 my-mut-effect-size
    ;parent2:
    let mut-rate-parent2 [my-mut-rate] of my-partner
    let MV-parent2 [my-mut-effect-size] of my-partner

    ;DEBUG:
    ;write "mut-rate-parent1: " write mut-rate-parent1 write " and mean: " print mean mut-rate-parent1
    ;write "MV-parent1: " write MV-parent1 write " and mean: " print mean MV-parent1
    ;write "mut-rate-parent2: " write mut-rate-parent2 write " and mean: " print mean mut-rate-parent2
    ;write "MV-parent2: " write MV-parent2 write " and mean: " print mean MV-parent2

    ;store the mean of the distribution of mutation effects for each parent
    ;parent1
    let mu-parent1 mu-dist-effect-size
    ;parent2
    let mu-parent2 [mu-dist-effect-size] of my-partner

    ;let mean-parental-mut-rate my-mut-rate
    ;let mean-parental-mut-effect-size my-mut-effect-size
    ;DEBUG:
    ;write "mean-parental-mut-rate: " print mean-parental-mut-rate
    ;write "mean-parental-mut-effect-size: " print mean-parental-mut-effect-size

    ; calculate the fecundity of the breeding pair pair-fecundity as in Björklund et al.
    ; Björklund(2009) assumed the fitness of the pair as the sum of the fitness
    ; values of the two parents (wsum) which is equivalent to the sum of their fecundities.
    ; Given that wsum is used, the population is saturated above the carrying capacity.
    let fecundity-of-pair ( fecundity + [ fecundity ] of my-partner )

    ;DEBUG:
    ;write "my-fecundity: " print fecundity
    ;write "fecundity-of-partner: " print [fecundity] of my-partner
    ;write "fecundity-of-pair: " print fecundity-of-pair

    ; set the genetic variance. Only for the calculation of environmental variance
    ; according to the heritability
    let genetic-variance-of-parents variance list (genetic-component)
                                                  ([genetic-component] of my-partner)
    let my-genetic-variance 0

    if-else ( how-genetic-variance = "parental-level")
    [
      set my-genetic-variance (genetic-variance-of-parents); parental-level
                                                           ; (Bjoerklund2009)
    ]
    [; else,
      if-else (how-genetic-variance = "parameter")
      [ set my-genetic-variance genetic-variance ] ; as in Reed et al. (2011):
        ;the population-level additive genetic variance is an input parameter (set by user)
      [; else
        if-else (how-genetic-variance = "population-level")
        [
        ; method as in Vincenzi & Piotti (2014): genetic variance is the total
        ; additive genetic variance for the trait at the population level
          set my-genetic-variance pop-gv
        ]
        ; else, catch exception
        [ error "undefined method for genetic variance" stop ]
        ]
    ]
    ; additive genetic variance of offspring equals half the additive genetic variance:
    set my-genetic-variance (1 / 2) * (my-genetic-variance)

    let number-of-offspring random-poisson fecundity-of-pair
    ; DEBUG
    ;write "nr. offspring: " print number-of-offspring

    if (number-of-offspring >= 1) ; reproduction occurs
    [
      hatch number-of-offspring
      [
        ; set stage as juvenile
        set stage "juvenile"
        set color blue
        ; set random orientation of the newborn and move to a green patch
        set heading random 360
        move-to one-of patches with [pcolor = green]


        ;*********************************************************************************
        ; INHERITANCE GENETIC COMPONENT OF PHENOTYPE
        ;*********************************************************************************
        ; inheritance with recombination and mutations, explicit genetics:

        ;1st mutation takes place during the creation of gametes inside each parent
        ; Each locus is checked for mutations according to the given mutation rate.
        ; If the model is set for evolving mutation rate, each individual will have its
        ; own value of mutation rate and effect size of mutation affecting the inheritance
        ; process
        ; The effect-size of the mutation is implemeted using a normal distribution
        ; N(0, variance), where the variance controls for the effect-size (Vincenzi 2014).
        ;
        ; A mutation can be explicitly simulated per locus or as a single mutation per
        ; strain with probability number-of-loci * mut-rate-per-locus, as in Vincenzi (2014)
        ; the effect-size of mutations can be modified in the future. For example,
        ; Vincenzi (2014) based the effect-size of mutations on a normal distribution
        ; N(0, effect-size). Through modifications of the variance of the distribution
        ; one can control the effect size of the mutation. Note that using the normal
        ; distribution, the expected value is the mean (in the example, 0 effect). This
        ; is important considering that perhaps the effect of mutations is unpredictable
        ; following an uniform distribution.

        ; THE MUTATION RATE SHOULD BE IMPLEMENTED PER LOCUS AND PER PARENT!
        ;DEBUG:
        ;write "parent1 strain1: " print my-loci1
        ;write "parent1 strain2: " print my-loci2

        ;parent1
        let l 0 ; counter
        ; set MR depending on selected method for MR:
        let m-rate 0
        if-else (method-for-MR = "decimal")
        [ set m-rate mean mut-rate-parent1]
        [
          if-else (method-for-MR = "exponent")
          [ set m-rate 10 ^ (-1 * (mean mut-rate-parent1)) ]
          [ error "unspecified method for MR" stop  ]
        ]

        let var mean MV-parent1
        ; copy strains info to avoid modifications of the original variable
        let strain1-parent1 my-loci1
        let strain2-parent1 my-loci2
        while [l < number-of-loci]
        [
          if (random-float 1 < m-rate) ; mutation occurs
          [ ; dna-strain1 = my-loci1 = strain1-parent1

            let dummy random-normal mu-parent1 sqrt var

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            if([optimum] of patch-here < phenotype) ; optimum at the left (optimum push for smaller trait value)
            [ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the probability of beneficial mutations
            ] ; else if ([optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain1-parent1
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain1-parent1 )
            set strain1-parent1 replace-item l strain1-parent1 dummy

            ; DEBUG
            ;write "value after: " print item l strain1-parent1
          ]
          if (random-float 1 < m-rate) ; mutation occurs
          [ ; dna-strain2 = my-loci2 = strain2-parent1

            let dummy random-normal mu-parent1 sqrt var

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            if([optimum] of patch-here < phenotype) ; optimum at the left (optimum push for smaller trait value)
            [ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the probability of beneficial mutations
            ] ; else if ([optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain2-parent1
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain2-parent1 )
            set strain2-parent1 replace-item l strain2-parent1 dummy

            ; DEBUG
            ;write "value after: " print item l strain2-parent1
          ]
          set l (l + 1)
        ]; while loop

        ;parent2
        set l 0 ; counter
        set m-rate 0
        if-else (method-for-MR = "decimal")
        [ set m-rate mean mut-rate-parent2]
        [
          if-else (method-for-MR = "exponent")
          [ set m-rate 10 ^ (-1 * (mean mut-rate-parent2)) ]
          [ error "unspecified method for MR" stop  ]
        ]

        set var mean MV-parent2
        let strain1-parent2 partner-loci1
        let strain2-parent2 partner-loci2
        while [l < number-of-loci]
        [
          if (random-float 1 < m-rate) ; mutation occurs
          [ ; dna-strain1 = partner-loci1= strain1-parent2

            let dummy random-normal mu-parent2 sqrt var

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            if([optimum] of patch-here < phenotype) ; optimum at the left (optimum push for smaller trait value)
            [ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the probability of beneficial mutations
            ] ; else if ([optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain1-parent2
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain1-parent2 )
            set strain1-parent2 replace-item l strain1-parent2 dummy

            ; DEBUG
            ;write "value after: " print item l strain1-parent2
          ]
          if (random-float 1 < m-rate) ; mutation occurs
          [ ; dna-strain2 = partner-loci2 = strain2-parent2

            let dummy random-normal mu-parent2 sqrt var

            ; a mutation in the direction of the optimum occurs with probability p =
            ; % beneficial mutations / 100
            if([optimum] of patch-here < phenotype) ; optimum at the left (optimum push for smaller trait value)
            [ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
                                     ; account for the probability of beneficial mutations
            ] ; else if ([optimum] of patch-here > phenotype) no need of further modification

            ; DEBUG
            ;write "value before: " print item l strain2-parent2
            ;write "effect of mut: " print dummy

            set dummy ( dummy + item l strain2-parent2 )
            set strain2-parent2 replace-item l strain2-parent2 dummy

            ; DEBUG
            ;write "value after: " print item l strain2-parent2
          ]
          set l (l + 1)
        ]; while loop

        ; 2nd inheritance with recombination:
        let recombining-allele-p1 0 ; 1 for strain1; 2 for strain2 of parent1
        let recombining-allele-p2 0 ; 1 for strain1; 2 for strain2 of parent2

        set l 0 ; counter
        set dna-strain1 n-values number-of-loci [0] ; initialization
        set dna-strain2 n-values number-of-loci [0]
        while [l < number-of-loci]
        [
          ; offspring strain 1, from parent 1
          set dna-strain1 replace-item l dna-strain1 (one-of (list item l strain1-parent1
                                                                   item l strain2-parent1))
          ; identify the recombining allele of parent1
          if-else (item l dna-strain1 = item l strain1-parent1)
          [ set recombining-allele-p1 1                     ] ; if strain1
          [ set recombining-allele-p1 2                     ] ; if strain2

          ; DEBUG:
          ;type item l dna-strain1 type " = " print item l strain1-parent1
          ;type item l dna-strain1 type " = " print item l strain2-parent1
          ;print recombining-allele-p1

          ; offspring strain 2, from parent 2
          set dna-strain2 replace-item l dna-strain2 (one-of (list item l strain1-parent2
                                                                   item l strain2-parent2))
          ; identify the recombining allele of parent2
          if-else (item l dna-strain2 = item l strain1-parent2)
          [ set recombining-allele-p2 1                     ] ; if strain1
          [ set recombining-allele-p2 2                     ] ; if strain2

          ; DEBUG:
          ;type item l dna-strain2 type " = " print item l strain1-parent2
          ;type item l dna-strain2 type " = " print item l strain2-parent2
          ;print recombining-allele-p2

          set l (l + 1)
        ] ; while loop
        ; DEBUG
        ;write "offspring strain1: " print dna-strain1
        ;write "offspring strain2: " print dna-strain2

        ; end of inheritance explicit genetics
        ;***********************************************************************************
        ;*********************************************************************************
        ; INHERITANCE MUTATION RATE AND MUTATION EFFECT SIZE
        ;*********************************************************************************
        ; inheritance of mutation rate and mutation effect size with probability of
        ; mutation
        ; DEBUG:
        ;print "before mutations"
        ;write "mutator-locus parent1: " print mut-rate-parent1
        ;write "effect-size parent1: " print MV-parent1
        ;write "mutator-locus parent2: " print mut-rate-parent2
        ;write "effect parent2: " print MV-parent2

        ; store original values
        let mut-backup-parent1 mut-rate-parent1
        let MV-backup-parent1 MV-parent1
        let mut-backup-parent2 mut-rate-parent2
        let MV-backup-parent2 MV-parent2
        let min-MR 0
        if(method-for-MR = "decimal")
        [ set min-MR 0.0001 ] ; minimum value of mutation rate 10^-4
        if(method-for-MR = "exponent")
        [ set min-MR 4 ]

        ; CHECK FOR MUTATIONS
        if (evolving-mut-rate? = "true" or
            evolving-mut-rate? = "only-mut-rate")
        [
          ; loop initialization

          let i 0
          while [ i < length mut-rate-parent1 ]
          [
            ;Parent1:
            ; in Cobben et al (2017) all loci, including the mutator locus are affected
            ; by the same mutation rate and effect size. Thus,
            set m-rate 0
            if-else (method-for-MR = "decimal")
            [ set m-rate mean mut-rate-parent1]
            [
              if-else (method-for-MR = "exponent")
              [ set m-rate 10 ^ (-1 * (mean mut-rate-parent1)) ]
              [ error "unspecified method for MR" stop  ]
            ]
            set var mean MV-parent1
            ;type "m-rate: " print m-rate
            let rn random-float 1
            ;if (rn < 0.0001)
            ;[type "rn: " print rn]

            ;if (random-float 1 < GA-m-rate) ; old method, unique mutation rate
            if (rn < m-rate) ; mutation occurs
            [
             ; print "mutation occurred"
              ;let dummy random-normal 0 sqrt GA-variance ; old method, unique effect size
              let dummy random-normal 0 sqrt var
              set dummy (dummy + item i mut-rate-parent1)
              set mut-rate-parent1 replace-item i mut-rate-parent1 dummy
              ; check for values out of range
              if(method-for-MR = "decimal")
              [
                if (item i mut-rate-parent1 > 1)
                [set mut-rate-parent1 replace-item i mut-rate-parent1 1]
                if( item i mut-rate-parent1 < min-MR)
                [set mut-rate-parent1 replace-item i mut-rate-parent1 min-MR] ; min mut rate 10^-4
              ]
              if(method-for-MR = "exponent")
              [ if( item i mut-rate-parent1 > min-MR)
                [set mut-rate-parent1 replace-item i mut-rate-parent1 min-MR] ; min mut rate 10^-4
                if( item i mut-rate-parent1 < 0 )
                [set mut-rate-parent1 replace-item i mut-rate-parent1 0]
              ]
            ]
            ;Parent2:
            ; in Cobben et al (2017) all loci, including the mutator locus are affected
            ; by the same mutation rate and effect size. Thus,

            set m-rate 0
            if-else (method-for-MR = "decimal")
            [ set m-rate mean mut-rate-parent2 ]
            [
              if-else (method-for-MR = "exponent")
              [ set m-rate 10 ^ (-1 * (mean mut-rate-parent2)) ]
              [ error "unspecified method for MR" stop   ]
            ]
            set var mean MV-parent2

            ;if (random-float 1 < GA-m-rate) ; old method, unique mutation rate
            if (random-float 1 < m-rate)
            [
              ;let dummy random-normal 0 sqrt GA-variance ;old method, unique effect size
              let dummy random-normal 0 sqrt var
              set dummy (dummy + item i mut-rate-parent2)
              set mut-rate-parent2 replace-item i mut-rate-parent2 dummy
              ; check for values out of range
              if(method-for-MR = "decimal")
              [
                if (item i mut-rate-parent2 > 1)
                [set mut-rate-parent2 replace-item i mut-rate-parent2 1]
                if( item i mut-rate-parent2 < min-MR)
                [set mut-rate-parent2 replace-item i mut-rate-parent2 min-MR] ; min mut rate 10^-4
              ]
              if(method-for-MR = "exponent")
              [ if( item i mut-rate-parent2 > min-MR)
                [set mut-rate-parent2 replace-item i mut-rate-parent2 min-MR] ; min mut rate 10^-4
                if( item i mut-rate-parent2 < 0 )
                [set mut-rate-parent2 replace-item i mut-rate-parent2 0]
              ]
            ]
            set i i + 1
          ]
        ]
        ; inheritance of mutation effect-size
        if (evolving-mut-rate? = "true" or
            evolving-mut-rate? = "only-mut-effect-size")
        [
          ; loop initialization
          let i 0
          while [ i < length MV-parent1]
          [
            ;Parent1:
            if (random-float 1 < GA-m-rate) ; mutation occurs
            [
              let dummy random-normal 0 sqrt GA-variance
              set dummy (dummy + item i MV-parent1)
              set MV-parent1 replace-item i MV-parent1 dummy
              ; check for values out of range
              if (item i MV-parent1 < 0)
              [set MV-parent1 replace-item i MV-parent1 0]
            ]
            ;Parent2:
            if (random-float 1 < GA-m-rate) ; mutation occurs
            [
              let dummy random-normal 0 sqrt GA-variance
              set dummy (dummy + item i MV-parent2)
              set MV-parent2 replace-item i MV-parent2 dummy
              ; check for values out of range
              if (item i MV-parent2 < 0)
              [set MV-parent2 replace-item i MV-parent2 0]
            ]
            set i i + 1
          ]
        ]
        ;DEBUG:
        ;print "after mutation:"
        ;write "mutator-locus parent1: " print mut-rate-parent1 ;write " and mean: "
        ;print mean mut-rate-parent1
        ;write "effect-size parent1: " write MV-parent1 write " and mean: "
        ;print mean MV-parent1
        ;write "mutator-locus-parent2 " print mut-rate-parent2
        ;write "effect parent2: " print MV-parent2

        ; RECOMBINATION
        ; with probability of recombination: 0 unlinked to 1 complete linkage

        ; initialization
        set my-mut-rate n-values length mut-rate-parent1 [0]
        set my-mut-effect-size n-values length MV-parent1 [0]
        ; inheritance of mutator locus
        let i 0 ; inheritance from parent1

        ; probability of recombination
        let q 0.5 + (0.5 * probability-recombination) ; probability of recombining for allele in strain1

        if (recombining-allele-p1 = 2) ; if strain2 of parent1 is recombining
        [ set q (1 - q)  ] ; probability of no recombining for allele in strain1 of parent1

        ; DEBUG
        ;type "q " print q

        if-else ( random-float 1 < q)
        [ set my-mut-rate replace-item i my-mut-rate item 0 mut-rate-parent1 ]
        [ set my-mut-rate replace-item i my-mut-rate item 1 mut-rate-parent1 ]
        ;set my-mut-rate replace-item i my-mut-rate (one-of mut-rate-parent1) ; old method

        ;DEBUG:
        ;type item i my-mut-rate type " = " print item 0 mut-rate-parent1
        ;type item i my-mut-rate type " = " print item 1 mut-rate-parent1

        set i i + 1; inheritance from parent2

        ; probability of recombination
        set q 0.5 + (0.5 * probability-recombination) ; probability of recombining for allele in strain1

        if (recombining-allele-p2 = 2) ; if strain2 of parent2 is recombining
        [ set q (1 - q)  ] ; probability of no recombining for allele in strain1 of parent2

        ; DEBUG
        ;type "q " print q

        if-else ( random-float 1 < q)
        [ set my-mut-rate replace-item i my-mut-rate item 0 mut-rate-parent2 ]
        [ set my-mut-rate replace-item i my-mut-rate item 1 mut-rate-parent2 ]

        ;DEBUG:
        ;type item i my-mut-rate type " = " print item 0 mut-rate-parent2
        ;type item i my-mut-rate type " = " print item 1 mut-rate-parent2
        ;set my-mut-rate replace-item i my-mut-rate (one-of mut-rate-parent2) ; old method

        ; inheritance of effect-size
        set i 0
        set my-mut-effect-size replace-item i my-mut-effect-size (one-of MV-parent1)
        set i i + 1
        set my-mut-effect-size replace-item i my-mut-effect-size (one-of MV-parent2)

        ; UNDO CHANGES TO MUT-RATE AND EFFECT-SIZE OF PARENTS
        ; for next offspring
        set mut-rate-parent1 mut-backup-parent1
        set MV-parent1 MV-backup-parent1
        set mut-rate-parent2 mut-backup-parent2
        set MV-parent2 MV-backup-parent2


        ;DEBUG:
        ;print "after recombination:"
        ;print "offspring:"
        ;write "mutator-locus parent1: " write my-mut-rate write " mean: " print mean my-mut-rate
        ;write "effect-size parent1: " write my-mut-effect-size write " mean: " print mean my-mut-effect-size
        ;*********************************************************************************

        ; update the genetic component of the offspring
        set genetic-component (sum dna-strain1) + (sum dna-strain2)

        ; DEBUG
        ;type "genetic-component: " print genetic-component

        ; environmental component 'e':
        ; Set the variance of environmental effect ve according to method of heritability
        ; if h2 = fixed, compute the variance of environmental effect according to the
        ; genetic variance and level of heritability
        ; else, ve is a constant parameter
        let ve 0 ; initialization
        if-else (how-heritability? = "fixed-value")
        [ set ve compute-env-effect-variance (heritability) (my-genetic-variance)]
        [ set ve env-effect-variance                                             ]

        set environmental-effect random-normal env-effect-mean sqrt ve

        ; DEBUG:
        ;write "genetic-variance: " print my-genetic-variance
        ;write "env-effect-variance: " print ve
      ] ; end of inheritance

   ] ; end of hatching n offspring ( reproduction)

    ; update reproduced? status of the partner to true
    ask my-partner [ set reproduced? true ]

  ]; end of if statement: whether a partner is available for reproduction

  ; update reproduced? of this turtle or individual to true:
    set reproduced? true

end





;**************************************************************************************
;        TO UPDATE-ENVIRONMENT
;**************************************************************************************
; Currently, there are two main environmental scenarios: climate-change and cyclic,
; respectively.
; In climate change scenario,
; the mean environmental optimum changes at rate r every iteration. For the model, each
; iteration is equivalent to a generation
; The optimum Q updates according to:
; 1) Qt = Q0 + rt (Directional deterministic)
; 2) Qt = Q0 + rt; Qt* = Qt + E, (Directional stochastic)
; where E is the type of noise defined as:
; Et = aE(t-1) + b*Dt, where D = N(0, 1)
; Here b is assumed to be b = sqrt VEt*(1 - a^2), as in Schwager et al (2006), where
; VEt is the environmental variance at time t; and the
; parameter a (set by the user), determines the strength of the autocorrelation
; (based on Ripa and Lundberg 1996, and Schwager et al. 2006)
;      a = 0 (white noise)
;  0 < a < 1 (red noise) Björklund et al used a = 0.7
; -1 < a < 0 (blue noise) Björklund et al used a = -0.7
;
; The parameter VEt can change or not in time, depending on the rate of change k of
; the variance. Thus, VEt = VE + kt, where VE is the inital environmental variance.
; This consideration allows to simulate the increased probability of extreme events
; as predicted by climatic IPCC scenarios, as in Vincenzi et al (2012).
; The function works in a patch context, example: ask patches [ update-optimum ]
;
; In the cyclic environment,
; The optimum Qt changes according to a sinusoidal function that considers two
; parameters:
; the amplitude A, and the period T, thus:
; Qt = A.sin(2.pi.t / T), according to Burger and Krall (2004)
to update-environment

  ; initial optimum of the environment Q0 = tita-0 according to Björklund2009
  let tita-0 0
  let new-optimum 0
  let new-noise 0
  let VE env-variance
  let k rate-change-of-env-variance
  let autocorr level-autocorr
  let b sqrt (1 - (autocorr ^ 2)) ; scaling factor of the variance according to
                           ; Schwager et al (2006)

  ; first catch exception:
  if (scenario? != "climate-change" and scenario? != "cyclic")
  [ error "update-environment: undefined scenario of environment" stop ]

  ; update the optimum according to the scenario of environment
  if (scenario? = "climate-change")
  [
    ; set parameter values for the scenario climate-change
    let r rate-change-of-optimum

    ; update moving mean envirnomental optimum
    if (ticks > time-to-balance)
    [ set new-optimum tita-0 + (r * (ticks - time-to-balance) ) ]; 1) directional deterministic
    if(ticks > (time-limit))
    [ set new-optimum  tita-0 + (r * (time-limit - time-to-balance) )] ; there is no longer a trend in mean environment
  ]

  if (scenario? = "cyclic") ; according to Burger and Krall 2004
  [
    ; set parameter values for the cyclic scenario
    let A 1 ; amplitude of the wave
    let T 1 ; period of the wave T = 1 / fr; fr is the frequency

    ; update the environmental optimum
    set new-optimum ( A * sin ( (2 * pi * (ticks - time-to-balance)) / T) )
  ]

  ; check for increasing environmetal variance and update the environmental variance
  if (ticks > time-to-balance)
  [set VE VE + (k * (ticks - time-to-balance) )]

  ; apply stochasticity around the environmental optimum
  ; stochastic environment with potentially different colour of noise
  ; depending on the level of outocorrelation a (level-autocorr)
  ;type "noise before " print noise
  set new-noise (autocorr * noise) + (b * (sqrt VE) * (random-normal 0 sqrt 1))
  set new-optimum new-optimum + new-noise
  ;type "noise after " print noise

  ; update optimum and noise of the patch
  set optimum new-optimum
  set noise new-noise

  ask patches with [distance myself <= size-of-environment]
  [
    set optimum new-optimum
    set noise new-noise
  ]

  ; DEBUG
  ;type "new optimum: " print optimum
  ;type "new noise: "   print noise
  ;type "patch optimum " print optimum
  ;type "patch noise: " print noise
  ;type "env-variance: " print env-variance
  ;type "new env-variance: " print VE

end






;**************************************************************************************
;        TO-REPORT COMPUTE-ENV-EFFECT-VARIANCE
;**************************************************************************************
; This function reports the corresponding value of environmental effect variance given the
; specified values of heritability and genetic variance
; h2: narrow-sense heritability
; gv: additive genetic variance
to-report compute-env-effect-variance [h2 gv]

  ; catch exception h2 <= 0:
  if-else (h2 <= 0)
  [ error "heritability muss be greater than 0 (h2 > 0)" ]
  [ report (gv / h2 ) - gv ] ; else

end






;**************************************************************************************
;        TO-REPORT COMPUTE-DEGREE-MALADAPTATION
;**************************************************************************************
; this function compute the degree of maladaptation according to Björklund et al (2009).
; Here this function is used to replicate the results in the mentioned work. They ran their
; simulation experiment for 30 generations and compute the evolutionary load, here degree
; of maladaptation, defined as:
; the sum(from i = 1 to t = 30) of [(average_z - opt)^2] / gamma
to-report compute-degree-maladaptation

 if-else ( count turtles > 0 )
 [
   let z mean [ phenotype ] of turtles
   let opt [ optimum ] of patch 0 0
   let gamma strength-selection

   ; DEBUG:
   ;write "mean z: " print z
   ;write "optimum: " print opt

   ; compute the degree-maladaptation
   report ( ( z - opt ) ^ 2 ) / gamma
 ]
 [ report 0 ] ; else (extinct population)

end





;**************************************************************************************
;        TO UPDATE-OUTPUT
;**************************************************************************************
; this function plots the distribution of phenotypes in the population
to update-output

  ;type "total " print count turtles
  ;type "adults " print count turtles with [stage = "adult"]
  ;type "juveniles " print count turtles with [stage = "juvenile"]

  ; plot time series:
  set-current-plot "Time series"
  set-current-plot-pen "turtles"
  plot count turtles

  ; plot genetic variance
  set-current-plot "genetic-variance"
  set-current-plot-pen "genetic-component"
  let dummy 0
  if(count turtles > 1)
  [ set dummy variance [genetic-component] of turtles ]
  plot dummy

  ; plot phenotypic variance in same plot of genetic-variance
  ;set-current-plot-pen "phenotype"
  ;set dummy 0
  ;if (count turtles > 1)
  ;[ set dummy variance [phenotype] of turtles ]
  ;plot dummy


  ; store the highest current phenotipic value in the population
  ; this value will be used below to set the range of the x axis
  ;let highest-value [phenotype] of max-one-of turtles [phenotype]
  set-current-plot "frequency distribution"


  ;set-plot-x-range -1 end-point-of-environment
  set-plot-pen-mode 1        ; (0 for lines, 1 for bars)
  set-histogram-num-bars 10
  ; plot the phenotypic frequency in the population:
  set-current-plot-pen "phenotype"
  histogram [phenotype] of turtles
  ; plot the genotypic frequency in the population:
  set-current-plot-pen "genotype"
  histogram [genetic-component] of turtles
  ; plot the environmental optimum
  set-current-plot-pen "env-optimum"
  plot-pen-reset
  set-plot-pen-mode 2
  ; plot a vertical line to show the current optimum phenotype in the environment
  let n 0
  while [n < carrying-capacity]
  [
    plotxy [optimum] of one-of patches with [pcolor = green] n
    set n n + 0.1
  ]

end






;**************************************************************************************
;        TO PLOT-MEAN-FITNESS
;**************************************************************************************
; plot mean fitness of turtles
to plot-mean-fitness

  ;type "total " print count turtles
  ;type "adults " print count turtles with [stage = "adult"]
  ;type "juveniles " print count turtles with [stage = "juvenile"]

  set-current-plot "degree of local adaptation"
  let dummy 0
  if (count turtles > 0 )
  [ set dummy mean [fitness] of turtles ]
  plot dummy

end






;**************************************************************************************
;        TO UPDATE-MEAN-DIST-EFFECTS
;**************************************************************************************
; update the mean of the distribution of effects according to the selected percentage
; of beneficial mutations and the variane of the distribution as given by the parameter
; "mutation-effect-size".
; this function operates in a turtle context
to update-mean-dist-effects

  ; compute the mean mu of the distribution of mutation effect size according to the
  ; proportions of beneficial mutations. Default values -> beneficial-mutations 0.5
  ; and mean 0. This mean "mu" is used to select the actual effect size during the
  ; process of inheritance after a mutation occurs
  set mu-dist-effect-size mean-dist (beneficial-mutations) (mean my-mut-effect-size)

  ;DEBUG:
  ; write "mu-dist-effects: " print mu-dist-effect-size

end






;**************************************************************************************
;        TO-REPORT MEAN-DIST
;**************************************************************************************
; this reporter returns the mean for the normal distribution "mean-dist"of mutation
; effect size according to the given proportion "p" of beneficial mutations (input
; parameter) and the mutation-effect-size or mutational variance "var". The mean-dist
; is used as the mean of the normal distribution from which a random variable is drawn.
; This reporter use the r-extension for Netlogo as it takes advantage from the qnorm
; function from the stats library in r
to-report mean-dist [p var]

  let x 0 ; point in the pdf. beneficial mutations are assumed to be any point => 0
  r:put "p" p
  r:put "sdev" sqrt var

  report x - (r:get "sdev * qnorm(1 - p, mean = 0, sd = 1)")

end






;**************************************************************************************
;        CHANGE LOG
;**************************************************************************************
;
;**************************************************************************************
;        MODEL VERSION: _V1.1
;**************************************************************************************
; the plot "mean fitness" is now called "degree of local adaptation"
; The functions implicit and explicit reproduction were modified in order to simulate
; lottery polygyny. This means that females reproduce only once, but males play a
; lottery
; The strength of selection of the fitness function were modified in order to yield
; similar results according to the following statements:
; a departure of one phenotypic deviation result in
; 60% of maximum fitness for a specialist (strongest strength of selection)
; 90% of max fitness for a moderate organism
; 95% of max fitness for a generalist (weakest level of strength of selection)
;**************************************************************************************
;        MODEL VERSION: _V1.2
;**************************************************************************************
; the mutation rate and mutation effect size can be set to initial random values (i.e.,
; individuals in the population differ initially in mutation rate and / or mutation
; effect size. There is an option to allow for fixed values as well
;**************************************************************************************
;        MODEL VERSION: _V1.2_experimental
;**************************************************************************************
; now stochasticity is simulated for the ciclic environment and if time-to-balance > 0
;**************************************************************************************
;        MODEL VERSION: _V2
;**************************************************************************************
; the user can now control the percentage of beneficial mutations (input parameter)
; a reporter that compute the corresponding mean for the normal distribution given the
; proportation of beneficial mutations "p" and the mutation effect size or mutational
; variance "var". This mean is then used to draw a random value from a normal
; distribution of mutation effect size
; The mean of the distribution of mutation effect size "mu" is
; computed at the beginning of the reproduction explicit genetics procedure using
; the mean-dist procedure
; The action of the mutator gene and its effect size is governed by the parental traits
; during the process of inheritance
;**************************************************************************************
;        MODEL VERSION: _V3
;**************************************************************************************
; the mean mu of the distribution of effect size is now a turtle own trait named
; mu-dist-effect-size.
; the update of this mu occurs through the function update-mean-dist-effects
; The process of inheritance in explicit genetics was changed such that, mutations
; take place before recombination. Mutations occur with probability as given by the
; mutator locus of the correspoding parental turtle
;**************************************************************************************
;        MODEL VERSION: _V3_1
;**************************************************************************************
; To Go procedure: update mean of the distribution of effect size only if % of
; beneficial mutation is not equal to 50% (i.e., 0.5), and if mutation-effect-size is
; evolving.
;  if ( beneficial-mutations != 0.5 )
;  [
;    if (evolving-mut-rate? = "true" or
;        evolving-mut-rate? = "only-mut-effect-size")
;    [ ask turtles [ update-mean-dist-effects ] ]
;  ]
;**************************************************************************************
;        MODEL VERSION: _V4
;**************************************************************************************
; Probability of recombination is now implemented. It works for any number of loci, and
; assumes that (when recombination probability > 0) the mutator locus is linked to the
; last locus affecting the ecological phenotype.
; All loci mutates with the same probability of mutation (as given by mutator locus),
; and same effect size, as in Cobben et al (2017)
;**************************************************************************************
;        MODEL VERSION: _V5
;**************************************************************************************
; The mutator locus can be simulated either as alleles of decimal values, or alleles
; coding for the exponent of base power 10 (as in Cobben et al 2017)
;**************************************************************************************
;        MODEL VERSION: _V33
;**************************************************************************************
; The distribution of mutation fitness effects according to different scenarios of
; beneficial mutations now works for all scenarios of environmental change. The following
; code was included:
; a mutation in the direction of the optimum occurs with probability p =
; % beneficial mutations / 100
;if([optimum] of patch-here < phenotype) ; optimum at the left (optimum push for smaller trait value)
;[ set dummy (-1) * dummy ; change sign such that the distribution of effects properly
;                         ; account for the probability of beneficial mutations
;] ; else if ([optimum] of patch-here > phenotype) no need of further modification
@#$#@#$#@
GRAPHICS-WINDOW
791
10
964
184
-1
-1
5.0
1
10
1
1
1
0
1
1
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
7
10
71
43
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
44
398
191
frequency distribution
NIL
freq
-10.0
10.0
0.0
10.0
true
true
"" ""
PENS
"phenotype" 1.0 0 -13791810 true "" ""
"env-optimum" 1.0 2 -2674135 true "" ""
"genotype" 1.0 0 -3026479 true "" ""

BUTTON
73
10
177
43
Go! one itera
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
179
10
295
43
Run the model!
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
191
398
338
Time series
time in generations
N
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"turtles" 1.0 0 -13791810 true "" ""

SLIDER
631
35
787
68
carrying-capacity
carrying-capacity
10
1000
1000.0
10
1
NIL
HORIZONTAL

SLIDER
631
69
787
102
population-size
population-size
10
1000
1000.0
10
1
NIL
HORIZONTAL

CHOOSER
631
103
787
148
type-organism
type-organism
"specialist" "moderate" "generalist"
1

CHOOSER
631
149
787
194
density-dependence-effect
density-dependence-effect
"weak" "normal" "strong" "very strong"
1

INPUTBOX
656
269
786
329
genetic-variance
0.2
1
0
Number

SLIDER
401
105
609
138
rate-change-of-optimum
rate-change-of-optimum
0
1
0.03
0.01
1
NIL
HORIZONTAL

CHOOSER
401
59
609
104
scenario?
scenario?
"climate-change" "cyclic"
0

SLIDER
512
315
654
348
heritability
heritability
0.1
1
1.0
0.1
1
NIL
HORIZONTAL

PLOT
7
338
398
486
genetic-variance
time in generations
genetic variation
0.0
2.0
0.0
1.0
true
false
"" ""
PENS
"genetic-component" 1.0 0 -16777216 true "" ""
"phenotype" 1.0 0 -3026479 true "" ""

CHOOSER
631
195
787
240
fitness-function
fitness-function
"Bjoerklund2009" "negative-exponential"
1

INPUTBOX
401
269
511
329
time-to-balance
100.0
1
0
Number

CHOOSER
656
329
786
374
how-genetic-variance
how-genetic-variance
"parameter" "parental-level" "population-level"
1

SLIDER
401
139
609
172
env-variance
env-variance
0
2
1.0
0.1
1
NIL
HORIZONTAL

PLOT
791
330
966
469
degree of local adaptation
time in generations
fitness
0.0
5.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

CHOOSER
401
330
511
375
how-genetics?
how-genetics?
"implicit" "explicit"
1

INPUTBOX
549
426
669
486
mut-rate-per-locus
0.0
1
0
Number

INPUTBOX
670
426
790
486
mut-effect-size
0.2
1
0
Number

SLIDER
549
394
790
427
number-of-loci
number-of-loci
1
50
1.0
1
1
NIL
HORIZONTAL

TEXTBOX
409
14
662
32
.................. ECOLOGY .................
14
54.0
1

TEXTBOX
401
247
785
281
..................................... EVOLUTION .....................................
14
93.0
1

SLIDER
297
10
398
43
time-limit
time-limit
10
1000
300.0
10
1
NIL
HORIZONTAL

SLIDER
401
207
609
240
level-autocorr
level-autocorr
-0.9
0.9
0.7
0.1
1
NIL
HORIZONTAL

CHOOSER
512
269
654
314
how-heritability?
how-heritability?
"fixed-value" "emerge"
0

SLIDER
401
173
609
206
rate-change-of-env-variance
rate-change-of-env-variance
0
0.1
0.0
0.01
1
NIL
HORIZONTAL

TEXTBOX
421
36
590
70
Environmental Scenario
14
54.0
1

TEXTBOX
643
16
793
34
Type of Organism
14
54.0
1

TEXTBOX
551
373
783
391
.............. Explicit genetics ..............
14
93.0
1

TEXTBOX
409
377
524
411
Implicit genetics
14
93.0
1

INPUTBOX
401
395
532
455
mutational-variance
0.01
1
0
Number

PLOT
791
187
964
327
environmental optimum
time
optimum
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -10899396 true "" "plot [optimum] of patch 0 0"

PLOT
969
209
1129
349
MR exponent
time
m-rate
0.0
1.0
0.0
0.01
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot 10 ^ (-1 * mean [mean my-mut-rate] of turtles)"

PLOT
969
349
1129
487
MR decimal
time
effect-size
0.0
1.0
0.0
0.01
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [mean my-mut-rate] of turtles"

CHOOSER
401
457
533
502
evolving-mut-rate?
evolving-mut-rate?
"false" "only-mut-rate" "only-mut-effect-size" "true"
1

INPUTBOX
972
146
1040
206
GA-m-rate
0.0
1
0
Number

INPUTBOX
1041
146
1126
206
GA-variance
0.0
1
0
Number

CHOOSER
972
54
1125
99
how-mut-rate?
how-mut-rate?
"random" "input-value"
0

CHOOSER
972
100
1125
145
how-mut-effect-size?
how-mut-effect-size?
"random" "input-value"
1

SLIDER
549
487
790
520
beneficial-mutations
beneficial-mutations
0.1
0.9
0.5
0.1
1
%
HORIZONTAL

SLIDER
791
471
965
504
probability-recombination
probability-recombination
0
1
1.0
0.1
1
NIL
HORIZONTAL

CHOOSER
972
10
1125
55
method-for-MR
method-for-MR
"exponent" "decimal"
0

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="this_is_it_r003" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean [MR] of turtles</metric>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;climate-change&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-heritability?">
      <value value="&quot;fixed-value&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-mut-rate?">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probability-recombination">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="-0.7"/>
      <value value="0"/>
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-m-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="method-for-MR">
      <value value="&quot;exponent&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beneficial-mutations">
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrying-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evolving-mut-rate?">
      <value value="&quot;only-mut-rate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-mut-effect-size?">
      <value value="&quot;input-value&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="300"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
