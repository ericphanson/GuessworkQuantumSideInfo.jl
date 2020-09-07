# A log of the output from the numeric experiments
```julia
GuessworkQuantumSideInfo.jl/numeric_experiments on  numeric_experiments [$?]
15:47:39 ❯ julia --project=.
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.5.1 (2020-08-25)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> versioninfo()
Julia Version 1.5.1
Commit 697e782ab8 (2020-08-25 20:08 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Core(TM) i7-6700K CPU @ 4.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-9.0.1 (ORCJIT, skylake)
Environment:
  JULIA_REVISE_POLL = 1

julia>

GuessworkQuantumSideInfo.jl/numeric_experiments on  numeric_experiments [$?] took 1m44s
15:49:29 ❯ julia --project=. --startup-file=no
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.5.1 (2020-08-25)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> include("run_problems.jl")
Academic license - for non-commercial use only
Problem 1/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:45
Finished!
[ Info: SubString{String}["1.254935008191061", "22.8479812"]
Problem 1/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:45
Finished!
[ Info: SubString{String}["1.254934996180558", "23.3493885"]
Problem 1/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:47
Finished!
[ Info: SubString{String}["1.2549349930201852", "24.7053996"]
Problem 1/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.254935008096305", "8.6902007"]
Problem 1/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.254934989188837", "8.458288"]
Problem 1/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4999999999980487", "8.7578906"]
Problem 1/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:32
Finished!
[ Info: SubString{String}["1.4999999999980487", "10.9102688"]
Problem 1/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4999999999980487", "9.4719107"]
Problem 1/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.254934948577658", "9.0020535"]
Problem 1/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.2549349952019324", "8.7551425"]
Problem 1/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.499999999879668", "9.0447568"]
Problem 1/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.499999999879668", "9.0737937"]
Problem 1/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.499999999879668", "9.0386222"]
Problem 2/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:56
Finished!
[ Info: SubString{String}["1.5315835908607798", "33.7600876"]
Problem 2/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:50
Finished!
[ Info: SubString{String}["1.5315835696345121", "27.0816645"]
Problem 2/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 2/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5315835982285744", "9.3482174"]
Problem 2/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5315835639649773", "9.1440335"]
Problem 2/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.55143345179782", "9.7268027"]
Problem 2/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.55143345179782", "9.6872898"]
Problem 2/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.55143345179782", "9.6855354"]
Problem 2/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5315836482762242", "9.0878709"]
Problem 2/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.5315835761159902", "8.8595117"]
Problem 2/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5514335194601956", "9.4359064"]
Problem 2/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5514335194601956", "9.4387474"]
Problem 2/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5514335194601956", "8.8855632"]
Problem 3/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:57
Finished!
┌ Warning: MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1529) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1530) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1531) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1532) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1533) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1534) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1535) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1536) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1537) of matrix 'A'.
│ MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1559) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1560) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1561) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1562) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1563) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1564) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1565) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1566) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1567) of matrix 'A'.
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:47
[ Info: SubString{String}["MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1529) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1530) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1531) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1532) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1533) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1534) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1535) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1536) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1537) of matrix 'A'.", "MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1559) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1560) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1561) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1562) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1563) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1564) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1565) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1566) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1567) of matrix 'A'.", "1.422649790197531", "32.9258842"]
Problem 3/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:50
Finished!
┌ Warning: MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1529) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1530) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1531) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1532) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1533) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1534) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1535) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1536) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1537) of matrix 'A'.
│ MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1559) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1560) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1561) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1562) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1563) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1564) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1565) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1566) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1567) of matrix 'A'.
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:47
[ Info: SubString{String}["MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1529) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1530) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1531) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1532) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1533) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1534) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1535) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1536) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1537) of matrix 'A'.", "MOSEK warning 705: #45 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1559) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1560) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1561) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1562) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1563) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1564) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1565) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1566) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(1567) of matrix 'A'.", "1.4226497382896444", "27.6751057"]
Problem 3/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 3/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.422649741094202", "9.326394"]
Problem 3/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.4226497338611328", "9.1115323"]
Problem 3/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.422649730078765", "9.7188591"]
Problem 3/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4226497077704807", "9.6625051"]
Problem 3/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4226497301069139", "9.8102977"]
Problem 3/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.4226494912208716", "8.976787"]
Problem 3/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.422649730963161", "8.741583"]
Problem 3/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4226497307011863", "9.2282523"]
Problem 3/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4226497307011854", "9.1781923"]
Problem 3/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4226497307011863", "9.3057953"]
Problem 4/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 4/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 4/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 4/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.709430585203046", "8.8681798"]
Problem 4/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.70943058402584", "8.6237742"]
Problem 4/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7094305753904342", "9.5134354"]
Problem 4/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7094305810180206", "9.7677162"]
Problem 4/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:32
Finished!
[ Info: SubString{String}["1.7094305248798227", "9.4909807"]
Problem 4/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7094306214229236", "8.5401641"]
Problem 4/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.7094305848837652", "8.2306411"]
Problem 4/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.7094305682408302", "9.0949739"]
Problem 4/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.7094305843004052", "9.2856356"]
Problem 4/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:28
Finished!
[ Info: SubString{String}["1.7094305682408737", "9.1506415"]
Problem 5/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:52
Finished!
[ Info: SubString{String}["1.255145077130398", "30.3665841"]
Problem 5/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:56
Finished!
[ Info: SubString{String}["1.255145070937713", "33.8677304"]
Problem 5/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 5/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.2551449846482652", "8.7806062"]
Problem 5/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.255144970810816", "8.5387358"]
Problem 5/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.4999999995692301", "8.8869137"]
Problem 5/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.4999999995692301", "8.8707072"]
Problem 5/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.4999999995692301", "8.9047032"]
Problem 5/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.2551450482354998", "8.4399617"]
Problem 5/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.2551449788825657", "8.3335215"]
Problem 5/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.499999997301028", "8.699171"]
Problem 5/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.499999997301028", "8.6579571"]
Problem 5/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.499999997301028", "8.7902263"]
Problem 6/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 6/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 6/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 6/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.471036475781681", "8.8248976"]
Problem 6/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:28
Finished!
[ Info: SubString{String}["1.4710363761428709", "8.5505169"]
Problem 6/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5562851836654978", "9.5011555"]
Problem 6/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.5562851836654978", "9.4520119"]
Problem 6/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.5562851836654978", "9.4534974"]
Problem 6/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:28
Finished!
[ Info: SubString{String}["1.4710366265761445", "8.4558049"]
Problem 6/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:28
Finished!
[ Info: SubString{String}["1.4710362238910173", "8.2458754"]
Problem 6/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5562846392846454", "9.6662718"]
Problem 6/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.5562846392846454", "9.107038"]
Problem 6/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.5562846392846454", "9.0146762"]
Problem 7/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 7/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 7/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 7/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:32
Finished!
[ Info: SubString{String}["1.7503204060173534", "8.8412297"]
Problem 7/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.750320341208066", "8.6292731"]
Problem 7/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.7763899365151006", "10.4054094"]
Problem 7/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7745147449590317", "10.4099024"]
Problem 7/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7763899365182758", "10.4019568"]
Problem 7/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7503204685792013", "8.8986134"]
Problem 7/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.7503217558950943", "8.3211818"]
Problem 7/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:32
Finished!
[ Info: SubString{String}["1.7763899165999122", "10.916788"]
Problem 7/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.7745140186833968", "10.2178801"]
Problem 7/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:31
Finished!
[ Info: SubString{String}["1.7744435625346329", "10.0230561"]
Problem 8/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:07
Timed out!
Problem 8/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:57
Finished!
[ Info: SubString{String}["2.7672182086517387", "32.2916971"]
Problem 8/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:45
Finished!
[ Info: SubString{String}["2.7676240967276766", "24.1154537"]
Problem 8/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:51
Finished!
[ Info: SubString{String}["2.7676240668569907", "30.4157592"]
Problem 8/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 8/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 8/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:58
Finished!
[ Info: SubString{String}["2.7672184622684357", "35.6556466"]
Problem 8/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:01:18
Finished!
[ Info: SubString{String}["2.7672181040034545", "54.6013249"]
Problem 8/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:01:37
Finished!
[ Info: SubString{String}["2.7672047842754983", "75.4263761"]
Problem 9/12, algorithm 1/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 2/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 3/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 4/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 5/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:09
Timed out!
Problem 9/12, algorithm 6/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:59
Finished!
[ Info: SubString{String}["2.320856634278685", "34.2465309"]
Problem 9/12, algorithm 7/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:53
Finished!
[ Info: SubString{String}["2.320856772167582", "30.996122"]
Problem 9/12, algorithm 8/13 100%|██████████████████████████████████████████████████████████████| Time: 0:00:51
Finished!
[ Info: SubString{String}["2.3208567801177336", "30.2602784"]
Problem 9/12, algorithm 9/13 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 10/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 9/12, algorithm 11/13 100%|█████████████████████████████████████████████████████████████| Time: 0:01:45
Finished!
[ Info: SubString{String}["2.320856591038295", "82.8359903"]
Problem 9/12, algorithm 12/13 100%|█████████████████████████████████████████████████████████████| Time: 0:01:13
Finished!
[ Info: SubString{String}["2.320856747573597", "50.7868997"]
Problem 9/12, algorithm 13/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:59
Finished!
[ Info: SubString{String}["2.3208566985412618", "37.362088"]
Problem 10/12, algorithm 1/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 10/12, algorithm 2/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 3/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 4/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:25
Finished!
┌ Error: ERROR: LoadError: OutOfMemoryError()
│ Stacktrace:
│  [1] Array at ./boot.jl:406 [inlined]
│  [2] Array at ./boot.jl:415 [inlined]
│  [3] similar at ./abstractarray.jl:675 [inlined]
│  [4] similar at ./abstractarray.jl:674 [inlined]
│  [5] _array_for at ./array.jl:678 [inlined]
│  [6] collect(::Base.Generator{UnitRange{Int64},GuessworkQuantumSideInfo.var"#10#15"{Int64}}) at ./array.jl:691
│  [7] guesswork(::Array{Float64,1}, ::Array{Array{Complex{Float64},2},1}; solver::MosekSolver, K::Int64, c::Array{Float64,1}, dual::Bool, remove_repetition::Bool, povm_outcomes::Combinatorics.MultiSetPermutations{Array{Int64,1}}, verbose::Bool, debug::Bool) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/src/SDP_formulation.jl:168
│  [8] macro expansion at ./timing.jl:310 [inlined]
│  [9] (::var"#23#28"{var"#21#26"})(::NamedTuple{(:p, :ρBs, :numeric_type, :problem),Tuple{Array{Float64,1},Array{Array{Complex{Float64},2},1},DataType,String}}) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/common.jl:79
│  [10] top-level scope at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
│  [11] include(::Function, ::Module, ::String) at ./Base.jl:380
│  [12] include(::Module, ::String) at ./Base.jl:368
│  [13] exec_options(::Base.JLOptions) at ./client.jl:296
│  [14] _start() at ./client.jl:506
│ in expression starting at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}[""]
ERROR: LoadError: BoundsError: attempt to access 1-element Array{SubString{String},1} at index [0:1]
Stacktrace:
 [1] throw_boundserror(::Array{SubString{String},1}, ::Tuple{UnitRange{Int64}}) at ./abstractarray.jl:541
 [2] checkbounds at ./abstractarray.jl:506 [inlined]
 [3] getindex at ./array.jl:815 [inlined]
 [4] run_problem(::Int64, ::Int64; verbose::Bool) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:51
 [5] top-level scope at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:61
 [6] include(::String) at ./client.jl:457
 [7] top-level scope at REPL[1]:1
in expression starting at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:60
```

Added error handling and continued...
```
julia> include("run_problems.jl")
Problem 10/12, algorithm 1/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 2/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 3/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 4/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:25
Finished!
┌ Error: ERROR: LoadError: OutOfMemoryError()
│ Stacktrace:
│  [1] Array at ./boot.jl:406 [inlined]
│  [2] Array at ./boot.jl:415 [inlined]
│  [3] similar at ./abstractarray.jl:675 [inlined]
│  [4] similar at ./abstractarray.jl:674 [inlined]
│  [5] _array_for at ./array.jl:678 [inlined]
│  [6] collect(::Base.Generator{UnitRange{Int64},GuessworkQuantumSideInfo.var"#10#15"{Int64}}) at ./array.jl:691
│  [7] guesswork(::Array{Float64,1}, ::Array{Array{Complex{Float64},2},1}; solver::MosekSolver, K::Int64, c::Array{Float64,1}, dual::Bool, remove_repetition::Bool, povm_outcomes::Combinatorics.MultiSetPermutations{Array{Int64,1}}, verbose::Bool, debug::Bool) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/src/SDP_formulation.jl:168
│  [8] macro expansion at ./timing.jl:310 [inlined]
│  [9] (::var"#23#28"{var"#21#26"})(::NamedTuple{(:p, :ρBs, :numeric_type, :problem),Tuple{Array{Float64,1},Array{Array{Complex{Float64},2},1},DataType,String}}) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/common.jl:79
│  [10] top-level scope at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
│  [11] include(::Function, ::Module, ::String) at ./Base.jl:380
│  [12] include(::Module, ::String) at ./Base.jl:368
│  [13] exec_options(::Base.JLOptions) at ./client.jl:296
│  [14] _start() at ./client.jl:506
│ in expression starting at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}[""]
┌ Error: Not enough returns, something went wrong
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:52
Problem 10/12, algorithm 5/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:24
Finished!
┌ Error: ERROR: LoadError: OutOfMemoryError()
│ Stacktrace:
│  [1] Array at ./boot.jl:406 [inlined]
│  [2] _array_for at ./array.jl:677 [inlined]
│  [3] guesswork(::Array{Float64,1}, ::Array{Array{Complex{Float64},2},1}; solver::MosekSolver, K::Int64, c::Array{Float64,1}, dual::Bool, remove_repetition::Bool, povm_outcomes::Combinatorics.MultiSetPermutations{Array{Int64,1}}, verbose::Bool, debug::Bool) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/src/SDP_formulation.jl:143
│  [4] macro expansion at ./timing.jl:310 [inlined]
│  [5] (::var"#24#29"{var"#21#26"})(::NamedTuple{(:p, :ρBs, :numeric_type, :problem),Tuple{Array{Float64,1},Array{Array{Complex{Float64},2},1},DataType,String}}) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/common.jl:82
│  [6] top-level scope at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
│  [7] include(::Function, ::Module, ::String) at ./Base.jl:380
│  [8] include(::Module, ::String) at ./Base.jl:368
│  [9] exec_options(::Base.JLOptions) at ./client.jl:296
│  [10] _start() at ./client.jl:506
│ in expression starting at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}[""]
┌ Error: Not enough returns, something went wrong
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:52
Problem 10/12, algorithm 6/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:47
Finished!
[ Info: SubString{String}["3.7811860979433263", "27.7171467"]
Problem 10/12, algorithm 7/13 100%|█████████████████████████████████████████████████████████████| Time: 0:01:29
Finished!
[ Info: SubString{String}["3.771192992150866", "66.469113"]
Problem 10/12, algorithm 8/13 100%|█████████████████████████████████████████████████████████████| Time: 0:02:35
Finished!
[ Info: SubString{String}["3.7682775167884057", "133.8114008"]
Problem 10/12, algorithm 9/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:24
Finished!
┌ Error: ERROR: LoadError: OutOfMemoryError()
│ Stacktrace:
│  [1] Array at ./boot.jl:406 [inlined]
│  [2] Array at ./boot.jl:415 [inlined]
│  [3] similar at ./abstractarray.jl:675 [inlined]
│  [4] similar at ./abstractarray.jl:674 [inlined]
│  [5] _array_for at ./array.jl:678 [inlined]
│  [6] collect(::Base.Generator{UnitRange{Int64},GuessworkQuantumSideInfo.var"#10#15"{Int64}}) at ./array.jl:691
│  [7] guesswork(::Array{Float64,1}, ::Array{Array{Complex{Float64},2},1}; solver::SCSSolver, K::Int64, c::Array{Float64,1}, dual::Bool, remove_repetition::Bool, povm_outcomes::Combinatorics.MultiSetPermutations{Array{Int64,1}}, verbose::Bool, debug::Bool) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/src/SDP_formulation.jl:168
│  [8] macro expansion at ./timing.jl:310 [inlined]
│  [9] (::var"#23#28"{var"#22#27"})(::NamedTuple{(:p, :ρBs, :numeric_type, :problem),Tuple{Array{Float64,1},Array{Array{Complex{Float64},2},1},DataType,String}}) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/common.jl:79
│  [10] top-level scope at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
│  [11] include(::Function, ::Module, ::String) at ./Base.jl:380
│  [12] include(::Module, ::String) at ./Base.jl:368
│  [13] exec_options(::Base.JLOptions) at ./client.jl:296
│  [14] _start() at ./client.jl:506
│ in expression starting at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}[""]
┌ Error: Not enough returns, something went wrong
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:52
Problem 10/12, algorithm 10/13 100%|████████████████████████████████████████████████████████████| Time: 0:00:24
Finished!
┌ Error: ERROR: LoadError: OutOfMemoryError()
│ Stacktrace:
│  [1] Array at ./boot.jl:406 [inlined]
│  [2] _array_for at ./array.jl:677 [inlined]
│  [3] guesswork(::Array{Float64,1}, ::Array{Array{Complex{Float64},2},1}; solver::SCSSolver, K::Int64, c::Array{Float64,1}, dual::Bool, remove_repetition::Bool, povm_outcomes::Combinatorics.MultiSetPermutations{Array{Int64,1}}, verbose::Bool, debug::Bool) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/src/SDP_formulation.jl:143
│  [4] macro expansion at ./timing.jl:310 [inlined]
│  [5] (::var"#24#29"{var"#22#27"})(::NamedTuple{(:p, :ρBs, :numeric_type, :problem),Tuple{Array{Float64,1},Array{Array{Complex{Float64},2},1},DataType,String}}) at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/common.jl:82
│  [6] top-level scope at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
│  [7] include(::Function, ::Module, ::String) at ./Base.jl:380
│  [8] include(::Module, ::String) at ./Base.jl:368
│  [9] exec_options(::Base.JLOptions) at ./client.jl:296
│  [10] _start() at ./client.jl:506
│ in expression starting at /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/do_problem.jl:20
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}[""]
┌ Error: Not enough returns, something went wrong
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:52
Problem 10/12, algorithm 11/13 100%|████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
[ Info: SubString{String}["3.779617362857275", "28.760714"]
Problem 10/12, algorithm 12/13 100%|████████████████████████████████████████████████████████████| Time: 0:01:32
Finished!
[ Info: SubString{String}["3.7715124469415007", "70.5699564"]
Problem 10/12, algorithm 13/13 100%|████████████████████████████████████████████████████████████| Time: 0:03:22
Finished!
[ Info: SubString{String}["3.768213428570261", "177.6444008"]
Problem 11/12, algorithm 1/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 2/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 3/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 4/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7368579207940624", "9.2295285"]
Problem 11/12, algorithm 5/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:29
Finished!
[ Info: SubString{String}["1.7368579132221822", "8.8348357"]
Problem 11/12, algorithm 6/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 7/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 8/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 9/13 100%|█████████████████████████████████████████████████████████████| Time: 0:00:32
Finished!
[ Info: SubString{String}["1.7368580657916646", "11.2161223"]
Problem 11/12, algorithm 10/13 100%|████████████████████████████████████████████████████████████| Time: 0:00:30
Finished!
[ Info: SubString{String}["1.7368575755119444", "9.1961901"]
Problem 11/12, algorithm 11/13 100%|████████████████████████████████████████████████████████████| Time: 0:02:02
Finished!
[ Info: SubString{String}["1.7554742953569498", "101.0598171"]
Problem 11/12, algorithm 12/13 100%|████████████████████████████████████████████████████████████| Time: 0:00:53
Finished!
[ Info: SubString{String}["1.7578147992900353", "32.9359894"]
Problem 11/12, algorithm 13/13 100%|████████████████████████████████████████████████████████████| Time: 0:00:44
Finished!
[ Info: SubString{String}["1.7583289930168684", "22.6860797"]
Problem 12/12, algorithm 1/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 12/12, algorithm 2/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 12/12, algorithm 3/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 12/12, algorithm 4/13 100%|█████████████████████████████████████████████████████████████| Time: 0:04:09
Finished!
[ Info: SubString{String}[""]
┌ Error: Not enough returns, something went wrong
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:52
Problem 12/12, algorithm 5/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 12/12, algorithm 6/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 12/12, algorithm 7/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 12/12, algorithm 8/13 100%|█████████████████████████████████████████████████████████████| Time: 0:05:02
Timed out!
Problem 12/12, algorithm 9/13 100%|█████████████████████████████████████████████████████████████| Time: 0:03:33
Finished!
[ Info: SubString{String}[""]
┌ Error: Not enough returns, something went wrong
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:52
Problem 12/12, algorithm 10/13 100%|████████████████████████████████████████████████████████████| Time: 0:05:03
Timed out!
Problem 12/12, algorithm 11/13 100%|████████████████████████████████████████████████████████████| Time: 0:01:47
Finished!
[ Info: SubString{String}["2.6455718400755233", "84.4330259"]
Problem 12/12, algorithm 12/13 100%|████████████████████████████████████████████████████████████| Time: 0:02:56
Finished!
[ Info: SubString{String}["2.639131319680283", "155.5140639"]
Problem 12/12, algorithm 13/13 100%|████████████████████████████████████████████████████████████| Time: 0:05:03
Timed out!
39×8 DataFrame. Omitted printing of 4 columns
│ Row │ algo                                │ settings                           │ numeric_type │ problem     │
│     │ String                              │ String                             │ String       │ String      │
├─────┼─────────────────────────────────────┼────────────────────────────────────┼──────────────┼─────────────┤
│ 1   │ MISDP                               │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ BB84(2)     │
│ 2   │ MISDP                               │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ BB84(2)     │
│ 3   │ MISDP                               │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ BB84(2)     │
│ 4   │ SDP                                 │ Mosek                              │ Float64      │ BB84(2)     │
│ 5   │ dual_SDP                            │ Mosek                              │ Float64      │ BB84(2)     │
│ 6   │ guesswork_upper_bound(max_time=20)  │ Mosek                              │ Float64      │ BB84(2)     │
│ 7   │ guesswork_upper_bound(max_time=60)  │ Mosek                              │ Float64      │ BB84(2)     │
│ 8   │ guesswork_upper_bound(max_time=240) │ Mosek                              │ Float64      │ BB84(2)     │
│ 9   │ SDP                                 │ SCS                                │ Float64      │ BB84(2)     │
│ 10  │ dual_SDP                            │ SCS                                │ Float64      │ BB84(2)     │
⋮
│ 29  │ MISDP                               │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ 3qutrits(2) │
│ 30  │ SDP                                 │ Mosek                              │ Float64      │ 3qutrits(2) │
│ 31  │ dual_SDP                            │ Mosek                              │ Float64      │ 3qutrits(2) │
│ 32  │ guesswork_upper_bound(max_time=20)  │ Mosek                              │ Float64      │ 3qutrits(2) │
│ 33  │ guesswork_upper_bound(max_time=60)  │ Mosek                              │ Float64      │ 3qutrits(2) │
│ 34  │ guesswork_upper_bound(max_time=240) │ Mosek                              │ Float64      │ 3qutrits(2) │
│ 35  │ SDP                                 │ SCS                                │ Float64      │ 3qutrits(2) │
│ 36  │ dual_SDP                            │ SCS                                │ Float64      │ 3qutrits(2) │
│ 37  │ guesswork_upper_bound(max_time=20)  │ SCS                                │ Float64      │ 3qutrits(2) │
│ 38  │ guesswork_upper_bound(max_time=60)  │ SCS                                │ Float64      │ 3qutrits(2) │
│ 39  │ guesswork_upper_bound(max_time=240) │ SCS                                │ Float64      │ 3qutrits(2) │
```

Added `MISDP_dB` setting and ran the problems with those options:
```
julia> include("run_problems.jl")
Problem 1/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:48
Finished!
[ Info: SubString{String}["1.2549350503546497", "22.9851179"]
Problem 1/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
[ Info: SubString{String}["1.2549350503546497", "23.6271632"]
Problem 1/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:48
Finished!
[ Info: SubString{String}["1.2549359025427322", "23.4738548"]
Problem 2/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:46
Finished!
[ Info: SubString{String}["1.531583578752242", "23.1706345"]
Problem 2/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:48
Finished!
[ Info: SubString{String}["1.5315835787522418", "23.630824"]
Problem 2/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:52
Finished!
[ Info: SubString{String}["1.5315876003496776", "27.1068016"]
Problem 3/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:48
Finished!
┌ Warning: MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(617) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(618) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(619) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(620) of matrix 'A'.
│ MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(629) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(630) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(631) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(632) of matrix 'A'.
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:47
[ Info: SubString{String}["MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(617) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(618) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(619) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(620) of matrix 'A'.", "MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(629) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(630) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(631) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(632) of matrix 'A'.", "1.4226497495723145", "23.2142919"]
Problem 3/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
┌ Warning: MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(617) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(618) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(619) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(620) of matrix 'A'.
│ MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(629) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(630) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(631) of matrix 'A'.
│ MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(632) of matrix 'A'.
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:47
[ Info: SubString{String}["MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(617) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(618) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(619) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(620) of matrix 'A'.", "MOSEK warning 705: #18 (nearly) zero elements are specified in sparse row ''(0) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(629) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(630) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(631) of matrix 'A'.", "MOSEK warning 705: #6 (nearly) zero elements are specified in sparse row ''(632) of matrix 'A'.", "1.4226497297957217", "23.5631458"]
Problem 3/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
┌ Error: ┌ Warning: Repeated integer solution without converging
│ └ @ Pajarito ~/.julia/packages/Pajarito/TFExZ/src/conic_algorithm.jl:1687
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}["1.422648765315336", "26.1466155"]
Problem 4/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
[ Info: SubString{String}["1.7094305955725198", "25.0106493"]
Problem 4/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:47
Finished!
[ Info: SubString{String}["1.7094305955725198", "23.5142924"]
Problem 4/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:01:39
Finished!
[ Info: SubString{String}["1.7094293534710712", "73.8001398"]
Problem 5/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:50
Finished!
[ Info: SubString{String}["1.255145017108577", "23.3100727"]
Problem 5/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:50
Finished!
[ Info: SubString{String}["1.255145017108577", "23.596305"]
Problem 5/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
[ Info: SubString{String}["1.2551468576010365", "24.7745932"]
Problem 6/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:49
Finished!
[ Info: SubString{String}["1.4710364763587662", "26.4991661"]
Problem 6/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:00:48
Finished!
[ Info: SubString{String}["1.4710364763587662", "24.4914356"]
Problem 6/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:01:57
Finished!
┌ Error: ┌ Warning: Repeated integer solution without converging
│ └ @ Pajarito ~/.julia/packages/Pajarito/TFExZ/src/conic_algorithm.jl:1687
└ @ Main /mnt/c/Users/eric/Code/GuessworkQuantumSideInfo.jl/numeric_experiments/run_problems.jl:42
[ Info: SubString{String}["1.4710321441071241", "95.0050658"]
Problem 7/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 7/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 7/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 8/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 1/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 3/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 9/12, algorithm 5/16 100%|██████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 1/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 3/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 10/12, algorithm 5/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 1/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 3/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 11/12, algorithm 5/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 12/12, algorithm 1/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 12/12, algorithm 3/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
Problem 12/12, algorithm 5/16 100%|█████████████████████████████████████████████████████████████| Time: 0:05:01
Timed out!
36×8 DataFrame. Omitted printing of 2 columns
│ Row │ algo     │ settings                           │ numeric_type │ problem     │ optval  │ elapsed_seconds │
│     │ String   │ String                             │ String       │ String      │ Float64 │ Float64         │
├─────┼──────────┼────────────────────────────────────┼──────────────┼─────────────┼─────────┼─────────────────┤
│ 1   │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ 2qubits(1)  │ 1.25494 │ 22.9851         │
│ 2   │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ 2qubits(1)  │ 1.25494 │ 23.6272         │
│ 3   │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ 2qubits(1)  │ 1.25494 │ 23.4739         │
│ 4   │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ 3qubits(1)  │ 1.53158 │ 23.1706         │
│ 5   │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ 3qubits(1)  │ 1.53158 │ 23.6308         │
│ 6   │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ 3qubits(1)  │ 1.53159 │ 27.1068         │
│ 7   │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ Y(1)        │ 1.42265 │ 23.2143         │
│ 8   │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ Y(1)        │ 1.42265 │ 23.5631         │
│ 9   │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ Y(1)        │ 1.42265 │ 26.1466         │
│ 10  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ BB84(1)     │ 1.70943 │ 25.0106         │
⋮
│ 26  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ Y(2)        │ NaN     │ NaN             │
│ 27  │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ Y(2)        │ NaN     │ NaN             │
│ 28  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ BB84(2)     │ NaN     │ NaN             │
│ 29  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ BB84(2)     │ NaN     │ NaN             │
│ 30  │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ BB84(2)     │ NaN     │ NaN             │
│ 31  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ 2qutrits(2) │ NaN     │ NaN             │
│ 32  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ 2qutrits(2) │ NaN     │ NaN             │
│ 33  │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ 2qutrits(2) │ NaN     │ NaN             │
│ 34  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=false) │ Float64      │ 3qutrits(2) │ NaN     │ NaN             │
│ 35  │ MISDP_dB │ Pajarito(Mosek, Gurobi, MSD=true)  │ Float64      │ 3qutrits(2) │ NaN     │ NaN             │
│ 36  │ MISDP_dB │ Pajarito(SCS, Cbc, MSD=false)      │ Float64      │ 3qutrits(2) │ NaN     │ NaN             │
```
