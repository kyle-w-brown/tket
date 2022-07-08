
Files in this "InitPlacement" directory are concerned with the Initial Qubit Placement problem.

Vaguely, this may be stated as:

INPUT:

    - A sequence of pairs of logical qubits, representing the 2-qubit gates in our abstract quantum circuit diagram, which we want to apply in that order.

    - A target graph: the architecture of a physical quantum computer: a weighted graph, where an edge means that a 2-qubit gate is possible between qubits situated at the endpoints, and the edge weight is proportional to the probability of gate error.

OUTPUT: an injective mapping f : {Logical qubits} -> {Physical qubits} (the initial qubit assignment).


NOTE ON GATE REORDERING: often the gates are reordered (e.g., partitioned into time slices, with disjoint gates occurring in parallel in each slice). But of course, even though they are parallel in the logical circuit, they might not be in the physical architecture, with token swapping taken into account. Thus there is still a lot of choice - even igoring the possibility of changing the GATES (with equivalent unitaries), there is still the problem of choosing the best REORDERING.

Of course, finding some f is trivial (since we haven't imposed any other constraints - we haven't yet got a subgraph monomorphism problem, weighted or unweighted).

What do we mean by a "good" f? Consider the following:

(1) It will, in general, be impossible to find f such that every desired 2-qubit gate can be achieved without extra qubit movement. So we MUST move qubits around in general, after the initial placement.

(2) We take the token swapping model: assume that extra swap gates will be added so as to move a pair of physical qubits into adjacent positions, when a gate between them needs to be applied.

(3) We assume that gate errors occur independently of each other (even on the same gate), and that the probability of error on a given qubit edge does not depend on WHICH gate is being applied. Thus we treat each swap gate as having the same error rate as any other 2-qubit gate on that edge (or maybe count a swap as 3 gates, if we treat it as 3 CNOTs; this is an adjustable parameter of course).

(4) We assume that single qubit gates have negligible error. (In future, single qubit error rates could also be taken into account, by modifying WSM to allow vertex weights as well as edge weights).

(5) We are trying to minimise the expected total number of 2-qubit gate errors (ideally, making it much less than 1!) This is simply (a constant times) the sum of the target edge weights over all gates applied (including any added swaps).

Now, there is vagueness in this; in (2), how will the swaps be calculated? There are many possibilities. Also, as stated above, two gates might be parallel in the logical circuit (acting on disjoint pattern edges), but added swap sequences to make them adjacent might have overlapping vertices, meaning that the gate order matters.

A proper solution would take all this into account, as a complete end-to-end discrete optimisation problem involving thousands of variables. (The token swaps included in the discrete variables, as well as all allowable gate reorderings). To be even fancier, you'd add in the possibilities of circuit identities to move and change gates.

However, this would obviously be extremely complex (I don't even know if there is any general-purpose algorithm or code in existence which could handle this well without extensive modification). One could probably encode everything as a massive integer programming problem, at the cost of introducing even more variables. This would also surely take a huge amount of computation even to find a good approximate solution.

So, our model with WSM is a HEURISTIC. The number of added swaps to make two vertices adjacent (ignoring other vertices) simply depends on the distance between them in the target graph.

However, we don't exactly do this to make the target weights; instead we combine target weights together recursively, and create a new target edge (u,v) whenever distance(u,v)=2, which then becomes like a new gate with higher error. We perform this through several generations.

Note that neither the pattern graph NOR the final target graph are inputs: the original weighted target graph (using gate errors) is expanded to a COMPLETE target graph by adding edges, and the pattern graph weights must be deduced from the gate sequence (which tells us not ONLY the edges - the interaction graph - but ALSO which edges are used more frequently than others, and so have a higher edge weight).

Logical qubit pairs with many gates are assigned HIGH edge weights, but desirable target edges (physical qubit pairs with low gate error rates) are assigned LOW edge weights.

The reason is, to minimise the sum of (pattern weight).(target weight) products, frequently used pattern edges will be driven towards lower target weights (we cannot avoid the pattern edges, but we can move them onto lower error gates, i.e. target edges with lower weights).

Also, the pattern weights have a decay factor (getting lower as the gates move further into the future - perhaps counting time using "time slices", rather than actual index in the sequence, but this is adjustable). (But note that we use arithmetic rather than geometric decay, since we want to avoid all non-integer arithmetic - we want identical results across different compilers and platforms, and floating point calculations are NOT necessarily identical across different platforms, compilers and "fast math" optimisation settings, etc.)

The reason is, we are doing only a rough HEURISTIC: as more swaps are added with some token swapping algorithm, more and more physical qubits will move further away from their original positions. Thus we want to concentrate more on the first few gates; it is less important to consider far future gates, as the qubits will be almost in "random" positions by the time we get to them.

Finally, we expect that other routines will exist, which are reasonable ITERATIVE algorithms: given an initial placement, try to improve it with something maybe similar to simulated annealing, with random jumping etc. For best results, our initial placement will only be the first input value to the iterative algorithm.

We stress once again that the above is all a HEURISTIC, it is NOT intended to be a particularly accurate model of how real token swapping etc. proceeds. But there are similar things in optimisation which work like this. E.g., in continuous optimisation, the BOBYQA algorithm works by constructing a quadratic model to the function being minimised. The model function is often quite a POOR approximation to the actual function; however, the SHAPE of the model function has lower values in roughly the same places, so that minimising it will still result in good values for the actual function.

We hope that the same will hold here; the rough heuristic function we are minimising will have a similar enough shape to the actual function to provide a good solution, EVEN THOUGH the actual detailed values will not be so good.

This is a rough first draft, there is a LOT of experimentation needed to find good parameter values. The difficult thing simulated annealing type algorithms is NOT the actual implementation, but getting something which deduces parameters fully automatically from the input problem data, which works well in most cases. Simulated annealing with human intervention, to solve a single given problem, is a lot easier.


MORE DETAILS - MCCT

We do something called MCCT (Monte Carlo Complete Target) to solve the initial WSM problem with COMPLETE target graph. This is like simulated annealing, but even simpler, since we only allow improving moves.

We just start from a random assignment (note that all assignments are POSSIBLE, just maybe bad, exactly because the target graph is complete. Because of this, WSM with complete target graph is far closer to the Travelling Salesperson problem than the more usual WSM problems).

We then do random jumps to improve the solution; if we don't decrease the scalar product fast enough, or find a solution (a "record breaker") better than the existing best known solution fast enough, we just switch to a new random assignment. (A better choice might be to perturb the best known solution randomly, but by a smaller amount, i.e. only a few random jumps/swaps. Something to try in future).

AFTER this solution is found, we then DELETE some unused target edges and vertices (actually, we start off with a new empty target, and add all vertices and edges used in the MCCT solution, plus some additional ones), to get a NEW WSM problem (but, since the new target graph is a subgraph of the complete graph, the MCCT solution is still a solution to the new WSM problem; and ANY solution to this new WSM problem also solves the WSM problem with complete target graph, with the same scalar product).

We then pass this new WSM problem into the existing WSM routines.


RESULTS SO FAR

We test with 2 fairly realistic models: binary trees, and square grids. So, we have edge weights representing gate error rates; we generate a random list of 2-qubit gates; and then, given any initial placement, we actually perform token swaps in a very simple way to enact the gates in the given order.

The point is NOT to benchmark tket's detailed existing routing; the point is to have a simple, but still reasonable, model.

What is the actual best possible placement, in this simplified model? We cannot compute it, obviously, but we just try running with many random placements, for several minutes, to see what kind of results are POSSIBLE, for comparison with our solutions generated by Initial Qubit Placement (IQP).

The results so far: for mapping, e.g. 30 or 90 qubits into a 10x10 square grid, or 30 qubits into a binary tree on 100 vertices, see the file:

tket/tests/WeightSubgrMono/InitPlacement/test_InitialPlacementProblems.cpp

It turns out in all these tests that the MCCT phase already finds a pretty good placement (in <50 ms), which takes ~1 minute or longer of random placements to beat.

This is a pretty good ratio, showing (at least) that the naive method of just trying random initial placements and augmenting with random jumping is likely to take a factor of 1000x or more to find an equally good placement.

Disappointingly, the main WSM routine does not improve on the MCCT solution at all in these tests (not because the WSM routines are bad, but because the MCCT seems to be ALREADY optimal for the new WSM problem produced).

However, we must stress that many more tests and experiments need to be done. The final result will probably be much better than this.


COMMENT ON PARTIAL WSM ASSIGNMENTS EXTENSION

One main idea was to use extension of partial solutions, i.e. add a callback to the main WSM solver, so that when the search node becomes nogood (further extension of the partial solution is impossible), the caller has the option of extending it to a full solution. This is because the internal main WSM-phase is solving a WSM problem with reduced target graph, but the earlier MCCT-phase involved a COMPLETE target graph. Thus, ANY partial assignment extends to a full assignment (as long as there is no target vertex clash).

[Notice how similar this is to graph colouring with unlimited colours, and Travelling Salesperson, rather than standard subgraph monomorphism. WSM seems to be quite powerful, in that it generalises BOTH standard subgraph monomorphism AND Travelling Salesperson, so I suspect it might be widely applicable, as a heuristic component for other problems].

To save time, we probably wouldn't extend EVERY partial solution, but maybe only some, according to some heuristic (such as, rough estimates for the final scalar product, given the number of PV, p-edges, and total p-weight so far).

Unfortunately, my crude extension algorithm seemed to be very poor (prioritising unassigned PV with high edge weight sum; assigning to TV with low edge weight sum). It seemed to be worse than just trying random assignments. Maybe it was wrong somehow, but it seems like we need a better extension algorithm to make this worthwhile. Maybe just crudely doing yet more random jumping is the way to go. (I.e., a cut-down version of MCCT on each extended partial solution).


COMMENT ON PARTIAL PLACEMENTS

Currently, if you only want to place some qubits, and other qubits have already been placed, there's no good way to tell the IQP routine that some target vertices are forbidden. The only way is to remove those vertices completely from the pattern and target graphs, but that's unsatisfactory because it also removes edges.

A better try is to modify MCCT and the main WSM routines to RESTRICT domains (i.e., specify at the start, for each PV, a SUBSET of TV within which Domain(PV) must lie).

This would actually be very easy to do, but I'm out of time. Definitely you should try this in future if you want to use IQP to place the PV in multiple chunks, rather than all at once.
