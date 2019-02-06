;;extensions [profiler]

globals
[   hours                     ;; records how many hours have passed after initiation of "go". 
    virusintegral             ;; approximation of integral of free virus
    infected-death-rate       ;; death rate of infected cell, corresponds to parameter d in paper
    max-time                  ;; maximum time to run simulation
    newvirions                ;; number of new virions an infected cell creates every hour
    total-targetcells         ;; total number of initial target cells 
    ctl-create-rate           ;; rate at which new CTL are created
    newctls                   ;; number of new ctls to potentially create
    firebreaks                ;; number of firebreaks (bystander cells killed by CTL)
    ctl-move-distance         ;; diffusion distance of CTL
    makemovie                  ;; set to zero or one depending on if we want to make a movie
]

breed [ctls a-ctl] ;; first breed is the lowest when drawing
breed [virions a-virion]
breed [uninfected-TC a-uninfected=TC]  ;; define different agents/breeds and give them both plural and singular names by which they are known to the program
breed [bystander-TC a-bystander-TC]
breed [infected-TC a-infected-TC]

ctls-own [ 
          boundstatus? ;;set to true if they are bound 
          ]

infected-TC-own
        [ infectage  ;;hours since infection
          boundstatus?
        ]

bystander-TC-own
        [ 
          boundstatus?
        ]

to setup ;; routine to setup simulation
 clear-all      ;; clear all
 setup-globals                  ;; initiates all global variables   
 setup-targetcells
 setup-infectedcells
 setup-virus
 setup-ctl
 if (makemovie = 1) 
     [
       movie-start "gjmovie.mov"
       movie-set-frame-rate 15
      ]
 ;;profiler:start         ;; start profiling
end  ;; end setup routine


to setup-globals
  set hours 0
  set newvirions 8 ;;number of virions created every hour
  set max-time 24 * 21 ;;maximum time of simulation
  set ctl-create-rate 1.0 ;;exponential growth rate of CTL
  set total-targetcells 10201
  set newctls 0
  set virusintegral 0
  set firebreaks 0
  set ctl-move-distance 2
  set makemovie 0
end

to setup-targetcells
  set-default-shape uninfected-TC "circle"
  set-default-shape bystander-TC "circle"
  set-default-shape infected-TC "circle"
  ask n-of total-targetcells patches
    [ sprout-uninfected-TC 1 
      [set color green]
    ]
end  

to setup-infectedcells
  ask uninfected-TC at-points [[0 0] ]
      [
        hatch-infected-TC 1  ;;create an infected cell
         [set color red
         set infectage 0
         set boundstatus? false]
        die ;;kill off uninfected cell
      ]
end  


to setup-virus
   create-virions 0
   [set color blue
    set heading random 360
   ]
end  

to setup-ctl
  set-default-shape ctls "dot"
  create-ctls 1 ;;initial CTL at time zero. 
   [set color white 
    setxy random-pxcor random-pycor
    set boundstatus? false
    ]
end  



to go ;; routine that runs the main simulation
  if ((count virions = 0) and (count infected-TC = 0))  
  [
   if (makemovie = 1) [movie-close]
   ;;set deadcells (total-targetcells - (count uninfected-TC))
   ;;profiler:stop          ;; stop profiling
   ;;print profiler:report  ;; view the results
   ;;profiler:reset         ;; clear the data
   stop
  ]
  ;;if (hours = max-time) [stop]  ;;could stop if max intregration time is reached
  if (makemovie = 1) [movie-grab-view]
  move-or-remove-virions  ;;move or remove free virions
  infect-cells-produce-bystander  ;;infect uninfected or bystander cells and/or create bystander cells
  form-complex-or-move-ctls ;;CTLs form complexes with infected/bystander cells or move
  kill-targetcells ;;CTL in complexes kill either infected or bystander cells
  make-new-virions-or-viruskill-cell ;; infected cells produce new virions and/or are killed by virus
  introduce-ctls ;;create new CTLs
  ;;if (hours = 150 and gj-rate = -5 and virion-move-distance = 0.2) [export-view "fig4a.png"]
  ;;if (hours = 150 and gj-rate = 0 and virion-move-distance = 0.2) [export-view "fig4b.png"]
  set virusintegral (virusintegral + (count virions) / 24)
  set hours (hours + 1) ;;timestep is 1h
end ;; end routine that runs main simulation

to move-or-remove-virions ;;first let virion move, then potentially kill it off
  ask virions
  [right random 360
  forward virion-move-distance
  if random 100 < virion-death-rate [die]  ;;any free virion can die with a certain probability 
  ]
end

to infect-cells-produce-bystander
  ask bystander-TC ;;first check all bystander cells to see if they might become infected
  [if any? virions-here
    [if random 100 < infect-rate
      [
        hatch-infected-TC 1  ;;create infected cell
         [set color red
         set infectage 0
         set boundstatus? false]
        die ;;kill off bystander cell
        ]
       ]
     ]
  ask uninfected-TC  ;;now check all uninfected cells. 
  [
    if any? infected-TC-on neighbors ;;first look if a neighboring infected cells might lead to bystander creation
      [if random-float 1 < 10 ^ gj-rate ;;gj-rate corresponds to parameter g in ODE model
        [
        hatch-bystander-TC 1
         [set color yellow
         set boundstatus? false]
         die 
        ]
       ]
    
    if any? virions-here   ;;now look if any virion here might lead to an infection. 
                           ;; note that bystander cells that are just created can get potentially immediately infected
     [if random 100 < infect-rate ;;infect-rate probability that a virion infects an uninfected or bystander cell, p_i
      [
        hatch-infected-TC 1  ;;create an infected cell
         [set color red
         set infectage 0
         set boundstatus? false]
        die ;;kill off uninfected cell
        ]
       ]
    ]
end
       
     


to introduce-ctls
  if (count ctls < 75000) ;;stop growing CTL beyond a certain number such that simulation can run in a reasonable time
   [set newctls  (newctls + ( ( ctl-create-rate / 24 ) * count ctls)) ;;approximate discrete version of exponential growth at rate ctl-create-rate
    if (newctls > 1)
      [create-ctls ceiling newctls 
        [set color white 
         setxy random-pxcor random-pycor
         set boundstatus? false
        ]
       set newctls 0
      ]
   ]
end

to form-complex-or-move-ctls
  ask ctls
  [
  if (count infected-TC-here > 0) 
   [if random 100 < complex-rate  ;; on-rate corresponds to p_i as defined in the paper
   [
    set boundstatus? true  ;;the CTL is tagged as bound
    ask infected-TC-here [set boundstatus? true]   ;;the present infected cell is tagged as bound
   ]
  ]
  if (count bystander-TC-here > 0) 
  [if random 100 < complex-rate 
   [set boundstatus? true
    ask bystander-TC-here [set boundstatus? true] 
   ]
  ]
 ;;if (count infected-TC-here = 0) and (count bystander-TC-here = 0)
   ;;  [set boundstatus? false]
    if not boundstatus?  ;;check if ctl are bound, only if not move them
    [right random 360
    forward ctl-move-distance
    ]
  ]
end


to kill-targetcells ;;killing occurs with some probability, depending on the time it takes to kill
  ask infected-TC
  [if (boundstatus? and random 100 < CTL-kill-rate)  ;;if the infected/bystander cells are bound by CTL, they will die. additionally, the CTL that were bound are released again.
    [
     ask CTLs-here [set boundstatus? false] 
     die
    ]
   ]
  ask bystander-TC
  [if (boundstatus? and random 100 < CTL-kill-rate) 
     [
      ask CTLs-here [set boundstatus? false] 
      set firebreaks (firebreaks + 1)
      die
     ]
  ]
end


to make-new-virions-or-viruskill-cell
  ask infected-TC
    [if infectage > -1 ;;wait at least this many ticks (hours, if we update every minute) before virus is produced. could be adjusted to introduce a latent phase
       [hatch-virions newvirions
        [set color blue
         set heading random 360
        ]
     ]
      ;;if random 100 < 4 [die] ;;die randomly at probability 0.04
      set infectage (infectage + 1)  ;;increase its age since infection by one. 
      if infectage = 24 [die] ;;die after 24 hours
    ]
end

@#$#@#$#@
GRAPHICS-WINDOW
303
10
1020
748
50
50
7.0
1
10
1
1
1
0
1
1
1
-50
50
-50
50
0
0
0
ticks

CC-WINDOW
5
762
1029
857
Command Center
0

BUTTON
163
476
226
509
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL

BUTTON
75
476
138
509
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL

MONITOR
105
372
166
417
virions
count virions
17
1
11

MONITOR
56
320
128
365
uninfected
count uninfected-TC
17
1
11

MONITOR
132
320
190
365
infected
count infected-TC
17
1
11

MONITOR
194
320
264
365
bystander
count bystander-TC
17
1
11

SLIDER
54
69
226
102
virion-death-rate
virion-death-rate
0
100
22
1
1
NIL
HORIZONTAL

SLIDER
55
147
227
180
gj-rate
gj-rate
-20
0
0
1
1
NIL
HORIZONTAL

MONITOR
46
372
103
417
ctls
count ctls
17
1
11

SLIDER
54
110
226
143
infect-rate
infect-rate
0
100
100
1
1
NIL
HORIZONTAL

SLIDER
55
188
227
221
complex-rate
complex-rate
0
100
100
1
1
NIL
HORIZONTAL

MONITOR
179
373
229
418
hours
hours
2
1
11

SLIDER
56
227
227
260
virion-move-distance
virion-move-distance
0
2
0.2
0.05
1
NIL
HORIZONTAL

SLIDER
52
26
224
59
CTL-kill-rate
CTL-kill-rate
0
100
100
1
1
NIL
HORIZONTAL

MONITOR
94
424
163
469
NIL
firebreaks
17
1
11

MONITOR
171
420
274
473
NIL
virusintegral
17
1
13

@#$#@#$#@
WHAT IS IT?
-----------
This model simulates an acute virus infection of a localized region, such as an area of epithelial tissue. The model focuses on gap-junction mediated antigen transport and CTL dynamics.

HOW IT WORKS
------------
For a description of the model, please see: "Gap junction mediated antigen transport and CTL responses during viral infections", Handel, Yates, Pilyugin and Antia. You can download this paper from the publications section on Andreas Handel's webpage: 
http://www.biology.emory.edu/Antia/ahandel/publications.htm 

The code itself (found under the 'Procedures' tab above) also contains comments that should make it easy to understand.


HOW TO USE IT
-------------
Setup and run the model. Sit back and watch :)


THINGS TO NOTICE
----------------
When the diffusion speed of the virions is slow, and the gap junction rate is high, it often occurs that the CTL kill all surrounding cells, leaving no cells the virus can enter in its vicinity. This is the firebreak mechanism described in the paper. If virions diffuse quickly, they can easily traverse any gap created by the CTL induced killing of target-cells, and the firebreak doesn't work anymore.


THINGS TO TRY
-------------
You can play with the different sliders to change the speed at which virus and CTL move, the rates of infection and complex formation, etc. If you want to change additional parameters, you can simply do so by editing the code, or introduce your own buttons/sliders.

EXTENDING THE MODEL
-------------------
Obviously, a viral infection is complicated, and many features of real infections have been neglected. One can easily extend this model to investigate other aspects of viral infections. For instance an interesting starting point might be to include some kind of chemical signal and allow CTL to follow this signal (chemotaxis) instead of diffusing around in a random fashion. Many more extensions are possible.


CREDITS AND REFERENCES
----------------------
Both the model and the accompanying paper can be found on Andreas Handel's webpage: 
http://www.biology.emory.edu/Antia/ahandel/index.htm

Feel free to use and modify the model as you wish. If you use the model or a modified version in your work, please cite this paper: "Gap junction mediated antigen transport and CTL responses during viral infections", Handel, Yates, Pilyugin and Antia (submitted)
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

link
true
0
Line -7500403 true 150 0 150 300

link direction
true
0
Line -7500403 true 150 150 30 225
Line -7500403 true 150 150 270 225

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

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="gj-sweep scenario three" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="63"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep slow diffusion" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep fast clearance" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep slow clearance" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep fast killing" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="86"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep slow killing" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="63"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very fast clearance-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep slow clearance-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="22"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep slow diffusion-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep fast diffusion-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep slow killing-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="63"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep fast killing-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="86"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very slow killing-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="39"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very fast killing-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep fast clearance-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very slow diffusion-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very slow clearance-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very fast diffusion-per fb" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-100"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="gj-sweep very fast clearance" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>virusintegral</metric>
    <metric>count uninfected-TC</metric>
    <metric>firebreaks</metric>
    <enumeratedValueSet variable="gj-rate">
      <value value="-4"/>
      <value value="-3"/>
      <value value="-2"/>
      <value value="-1"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infect-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="complex-rate">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-death-rate">
      <value value="34"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virion-move-distance">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CTL-kill-rate">
      <value value="100"/>
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
