import bisect
import math
import numpy as np

from .configuration import WormConfiguration
from .configuration import TYPE_HOP, TYPE_WORM_TAIL, TYPE_WORM_HEAD, TYPE_WORM_DUMMY, EPSILON

class WormAlgorithm:

    def __init__(
            self,
            lattice,
            hamiltonian,
            beta,
            n_max = 5,
            c_worm = 2.0,
            energy_off = 1.0,
            seed = None,
            epsilon_time = 1.0e-3,
            kink_prob = 0.4,
            glue_prob = 0.3,
            move_prob = 0.3,
            verbose = False,
            ):

        self.lattice = lattice
        self.hamiltonian = hamiltonian
        self.beta = beta
        self.n_max = n_max
        self.c_worm = c_worm
        self.energy_off = energy_off

        self.rng = np.random.default_rng(seed)
        self.config = WormConfiguration(lattice, hamiltonian, beta)

        self.stats = {
                'insert_attempts': 0,
                'insert_accepts': 0,
                'glue_attempts'  : 0,
                'glue_accepts'  : 0,
                'move_attempts' : 0,
                'move_accepts'  : 0,
                'insertkink_attempts' : 0,
                'insertkink_accepts'  : 0,
                'deletekink_attempts' : 0,
                'deletekink_accepts'  : 0,
                'sweeps'        : 0,
                }

        self.epsilon_time = epsilon_time * self.beta

        self.glue_prob = glue_prob
        self.move_prob = move_prob
        self.kink_prob = kink_prob
        self.verbose = verbose
        # helpers

    def _log(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)

    def _find_event_index(self, site, time, eps = EPSILON, type_filter = None):
        """
        Encontrar el indice siguiente que coincide con time dentro de la tolerancia
        """
        tolerance = eps * self.beta
        events = self.config.events[site]
        if not events:
            return None


        tnorm = self.config._norm_time(time)
        times = [e['time'] for e in events]

        i = bisect.bisect_left(times, tnorm)

        for index in (i - 1, i, i + 1):
            if 0 <= index < len(events):
                event = events[index]
                if self.config._time_close(event['time'], tnorm):
                    if type_filter is None or event['type'] == type_filter:
                        return index

        return None

    def _matrix_prod_from_occ_change(self, occ_before, occ_after):
        """
        matrix elements product <occ_before|b_j|occ_after><occ_after| b_j^dagger | occ_before>
        """
        b = self.hamiltonian.bosonic_matrix_element(occ_before, occ_after)
        b_dagger = self.hamiltonian.bosonic_matrix_element(occ_after, occ_before)
        return b * b_dagger

# Insert worm

    def insert_worm(self):
        """Switching between the partition function sector and the G func sector by inserting or removing a worm pair (head + tail).

        Can only be called form the z-sector"""

        if not self.config.in_z_sector:
            return False

        self.stats['insert_attempts'] += 1

        site = self.rng.integers(0, self.lattice.get_nsites()) #we choose a random site A
        time = self.rng.uniform(0.0, self.beta)                #At a random time \tau_A
        occ = self.config.get_occupation_at_time(site, time)   #We find the occupation at A
        new_occ = occ                                          #We add change the occupation later

        create_first = self.rng.random() < 0.5                 #We decide with equal probability whether we create of annihilate first
        if create_first:
            new_occ += 1
            if new_occ > self.n_max:
                return False                                   #Occupation out of boundaries (superior bound surpassed), we reject
        else :
            if occ == 0:                                        #Occupation out of boundaries (we can't annhilate particles if they do not exist), we reject
                return False

            new_occ -= 1

        mat_prod = self._matrix_prod_from_occ_change(occ, new_occ) #|< occ | b_A | new_occ >|^2

        if mat_prod <= 0.0:                                    #It shouldn't be less equal zero, but just in case
            return False

        acceptance_ratio = 2 * self.c_worm  * mat_prod         #Acceptance ratio
        acceptance_prob = min(1.0, acceptance_ratio)           #Acceptance probability
        if acceptance_prob < self.rng.random():                #if the probability is less than a random number, we reject
            return False

        epsilon = max(EPSILON * self.beta, EPSILON )           #We insert the worms a infinitesimal time appart
        time_tail = time
        time_head = time + epsilon

        index_tail = self.config.insert_element(site, time_tail, TYPE_WORM_TAIL, occ, new_occ, linked_site=site)
        index_head = self.config.insert_element(site, time_head, TYPE_WORM_HEAD, new_occ, occ, linked_site=site)

        self.config.worm_tail_site = site
        self.config.worm_tail_time = self.config.events[site][index_tail]['time']
        self.config.worm_head_site = site
        self.config.worm_head_time = self.config.events[site][index_head]['time']
        self.config.worm_head_wpm = 1
        self.config.worm_tail_wpm = -1
        self.config.in_z_sector = False
        self.stats['insert_accepts'] += 1
        self._log(f'[Insert Worm] site = {site}; times = ({time_tail:.6f}, {time_head:.6f});  occ:{occ} --> {new_occ}')
        return True

    def glue_worm(self, eps = EPSILON):
        if self.config.in_z_sector:
            return False

        self.stats['glue_attempts'] += 1

        head_site = self.config.worm_head_site
        tail_site = self.config.worm_tail_site

        # Check if head and tail are on the same site
        if head_site != tail_site:
            return False

        site = head_site
        time_head = self.config._norm_time(self.config.worm_head_time)
        time_tail = self.config._norm_time(self.config.worm_tail_time)

        # Find indices of head and tail events
        index_head = self._find_event_index(site, time_head, type_filter=TYPE_WORM_HEAD)
        index_tail = self._find_event_index(site, time_tail, type_filter=TYPE_WORM_TAIL)

        if index_head is None or index_tail is None:
            return False

        # Get events directly
        head_event = self.config.events[site][index_head]
        tail_event = self.config.events[site][index_tail]

        # Calculate matrix element product using events
        # Head event: occ_left -> occ_right (reverts tail change)
        # Tail event: occ_left -> occ_right (initial change)
        occ_before = tail_event['occ_left']
        occ_after = tail_event['occ_right']

        mat_prod = self._matrix_prod_from_occ_change(occ_before, occ_after)
        if mat_prod <= 0.0:
            return False

        acceptance_ratio = 1.0 / (2.0 * self.c_worm * mat_prod)
        acceptance_prob = min(1.0, acceptance_ratio)

        if acceptance_prob < self.rng.random():
            return False

        # Remove events (in reverse order to preserve indices)
        for index in sorted((index_head, index_tail), reverse=True):
            self.config.remove_element(site, index)

        # Reset worm state
        self.config.worm_head_site = -1
        self.config.worm_head_time = -1.0
        self.config.worm_tail_site = -1
        self.config.worm_tail_time = -1.0
        self.config.worm_head_wpm = 0
        self.config.worm_tail_wpm = 0
        self.config.in_z_sector = True

        self.stats['glue_accepts'] += 1
        self._log(f'[Glue] site = {site}; occ: {occ_before} --> {occ_after}')
        return True

    def move_worm(self, eps = EPSILON):

        tolerance = eps * self.beta
        if self.config.in_z_sector:           #Must be in the G sector
            return False

        self.stats['move_attempts'] += 1

        move_head = self.rng.random() < 0.5        #we decide with equal probability which worm we'll move
        forward = self.rng.random() < 0.5
        site, current_time, current_wpm, event_type = 0, 0.0, 0, 0 #initialization of values
        if move_head:
            site = self.config.worm_head_site
            current_time = self.config.worm_head_time
            current_wpm = self.config.worm_head_wpm
            event_type = TYPE_WORM_HEAD
        else:
            site = self.config.worm_tail_site
            current_time = self.config.worm_tail_time
            current_wpm = self.config.worm_tail_wpm
            event_type = TYPE_WORM_TAIL

        if current_wpm != 0:
            if (forward and current_wpm == 1) or (not forward and current_wpm == -1):
                return False #We reject movements towards existing elements

        E_L = self.config.compute_local_energy(site, current_time - eps)
        E_R = self.config.compute_local_energy(site, current_time + eps)


        rate, r_denom, r_num, pexp, proposed_time = 0.0, 0.0, 0.0, 0.0, 0.0 #Initialization

        if forward:
            tau = self.config.find_time_to_next_element(site, current_time, include_neighbors = True)

            if E_R > E_L:
                rate = self.energy_off
                pexp = -math.log(self.rng.random()) / rate
                r_denom = rate if tau > pexp else 1.0
                r_num = E_R - E_L + self.energy_off if current_wpm == 0 else 1.0
            else:
                rate = E_L - E_R + self.energy_off
                pexp = -math.log(self.rng.random()) / rate
                r_denom = rate if tau > pexp else 1.0
                r_num = self.energy_off if current_wpm ==  0 else 1.0
            proposed_time = current_time + min(pexp, tau)
        else:
            tau = self.config.find_time_to_prev_element(site,current_time, include_neighbors = True)
            if E_R > E_L:
                rate = E_R + E_L - self.energy_off
                pexp = -math.log(self.rng.random()) / rate
                r_denom = rate if tau > pexp else 1.0
                r_num = self.energy_off if current_wpm == 0 else 1.0
            else:
                rate = self.energy_off
                pexp = -math.log(self.rng.random()) / rate
                r_denom = rate if tau > pexp else 1.0
                r_num = E_L - E_R + self.energy_off if current_wpm == 0 else 1.0
            proposed_time = current_time - min(pexp, tau)

        proposed_time = self.config._norm_time(proposed_time)
        #Avoids roundoff errors in exchange of a very small discretization error.
        if pexp < tolerance or abs(tau - pexp) < tolerance:
            return False

        acceptance_ratio = r_num / r_denom
        acceptance_prob = min(1.0, acceptance_ratio)

        if (acceptance_prob < self.rng.random() ):
            return False

        index = self._find_event_index(site, current_time, type_filter = event_type)

        if index is None:
            return False

        event = self.config.events[site][index].copy()
        self.config.remove_element(site, index)
        new_index = self.config.insert_element(site, proposed_time, event['type'], event['occ_left'], event['occ_right'], event['linked_site'])

        new_wpm = 0 if pexp < tau else (-1 if forward else 1)
        if move_head:
            self.config.worm_head_wpm = new_wpm
            self.config.worm_head_time = proposed_time
        else:
            self.config.worm_tail_wpm = new_wpm
            self.config.worm_tail_time = proposed_time

        self.stats['move_accepts'] += 1
        self._log(f"[Move] site = {site}; type = {'H' if move_head else 'T' }; t_old = {current_time:.6f}; t_new={proposed_time:.6f}")
        return True


    def insert_kink(self):
        if self.config.in_z_sector:
            return False

        self.stats['insertkink_attempts'] += 1

        # Choose which worm end to move
        if self.rng.random() < 0.5:
            site = self.config.worm_head_site
            current_time = self.config.worm_head_time
            worm_type = TYPE_WORM_HEAD
            is_head = True
        else:
            site = self.config.worm_tail_site
            current_time = self.config.worm_tail_time
            worm_type = TYPE_WORM_TAIL
            is_head = False

        # Get neighbors
        neighbors = self.lattice.get_neighbors(site)
        if not neighbors:
            return False

        # Choose random neighbor
        neighbor_site = self.rng.choice(neighbors)

        # Get occupations at current time
        occ_site = self.config.get_occupation_at_time(site, current_time)
        occ_neighbor = self.config.get_occupation_at_time(neighbor_site, current_time)

        # Determine the transition
        if is_head:
            # Head represents annihilation: destroy at site, create at neighbor
            if occ_site < 1:
                return False  # Cannot annihilate if no particle
            new_occ_site = occ_site - 1
            new_occ_neighbor = occ_neighbor + 1
            mat_site = math.sqrt(occ_site)  # b|n> = √n|n-1>
            mat_neighbor = math.sqrt(occ_neighbor + 1)  # b†|n> = √(n+1)|n+1>
        else:
            # Tail represents creation: create at site, destroy at neighbor
            if occ_neighbor < 1:
                return False  # Cannot annihilate if no particle at neighbor
            new_occ_site = occ_site + 1
            new_occ_neighbor = occ_neighbor - 1
            mat_site = math.sqrt(occ_site + 1)  # b†|n> = √(n+1)|n+1>
            mat_neighbor = math.sqrt(occ_neighbor)  # b|n> = √n|n-1>

        # Check bounds
        if (new_occ_site < 0 or new_occ_site > self.n_max or
            new_occ_neighbor < 0 or new_occ_neighbor > self.n_max):
            return False

        # Acceptance ratio
        matrix_element = mat_site * mat_neighbor
        z = len(neighbors)
        # Probability to choose this specific move: 1/(2*z)
        # (1/2 for head/tail, 1/z for neighbor choice)
        proposal_prob_forward = 1.0 / (2.0 * z)

        # For reverse move (delete_kink), probability is 1 (deterministic if we're at the kink)
        proposal_prob_reverse = 1.0

        proposal_ratio = proposal_prob_reverse / proposal_prob_forward

        # Acceptance probability
        acceptance_ratio = (matrix_element * self.hamiltonian.t * self.beta *
                       proposal_ratio)
        acceptance_prob = min(1.0, acceptance_ratio)

        if acceptance_prob < self.rng.random():
            return False

        # Remove the old worm event
        worm_index = self._find_event_index(site, current_time, type_filter=worm_type)
        if worm_index is None:
            return False

        self.config.remove_element(site, worm_index)

        # Insert hopping events at the same time
        self.config.insert_element(site, current_time, TYPE_HOP,
                              occ_site, new_occ_site,
                              linked_site=neighbor_site)
        self.config.insert_element(neighbor_site, current_time, TYPE_HOP,
                              occ_neighbor, new_occ_neighbor,
                              linked_site=site)

        # Insert worm at neighbor site at SAME time (not shifted)
        # This is key: the worm should be at the same imaginary time as the hop
        index_worm = self.config.insert_element(neighbor_site, current_time,
                                           worm_type, new_occ_neighbor, occ_neighbor,
                                           linked_site=neighbor_site)

        # Update worm position
        if is_head:
            self.config.worm_head_site = neighbor_site
            self.config.worm_head_time = current_time
            self.config.worm_head_wpm = 0  # Reset wpm since worm is now at hop event
        else:
            self.config.worm_tail_site = neighbor_site
            self.config.worm_tail_time = current_time
            self.config.worm_tail_wpm = 0

        self.stats['insertkink_accepts'] += 1
        return True

    def delete_kink(self):
        if self.config.in_z_sector:
            return False

        self.stats['deletekink_attempts'] += 1
        site, current_time, current_wpm, worm_type = 0, 0.0, 0, 0

        move_head = self.rng.random() < 0.5

        if move_head:
            site = self.config.worm_head_site
            current_time = self.config.worm_head_time
            current_wpm = self.config.worm_head_wpm
            worm_type = TYPE_WORM_HEAD
        else:
            site = self.config.worm_tail_site
            current_time = self.config.worm_tail_time
            current_wpm = self.config.worm_tail_wpm
            worm_type = TYPE_WORM_TAIL

        # Allow wpm = 0, ±1 for kink deletion
        #if current_wpm != 0:
        #    return False

        # Find the hop event at current site
        hop_index = self._find_event_index(site, current_time, type_filter=TYPE_HOP)
        if hop_index is None:
            # Try with larger tolerance
            hop_index = self._find_event_index(site, current_time,
                                          eps=10*EPSILON, type_filter=TYPE_HOP)
            if hop_index is None:
                return False

        hop_event = self.config.events[site][hop_index]
        neighbor_site = hop_event['linked_site']

        # Find the corresponding hop event at neighbor site
        neighbor_hop_index = self._find_event_index(neighbor_site, current_time,
                                               type_filter=TYPE_HOP)
        if neighbor_hop_index is None:
            neighbor_hop_index = self._find_event_index(neighbor_site, current_time,
                                                   eps=10*EPSILON, type_filter=TYPE_HOP)
            if neighbor_hop_index is None:
                return False

        neighbor_hop_event = self.config.events[neighbor_site][neighbor_hop_index]

        if neighbor_hop_event['linked_site'] != site:
            return False

        occ_A = hop_event['occ_left']
        final_occ_A = hop_event['occ_right']
        occ_B = neighbor_hop_event['occ_left']
        final_occ_B = neighbor_hop_event['occ_right']

        # Calculate matrix elements based on the actual transition
        if final_occ_A == occ_A - 1 and final_occ_B == occ_B + 1:
            # This was a head-type transition: destroy at A, create at B
            mat_A = math.sqrt(occ_A)  # b|n> = √n|n-1>
            mat_B = math.sqrt(occ_B + 1)  # b†|n> = √(n+1)|n+1>
        elif final_occ_A == occ_A + 1 and final_occ_B == occ_B - 1:
            # This was a tail-type transition: create at A, destroy at B
            mat_A = math.sqrt(occ_A + 1)  # b†|n> = √(n+1)|n+1>
            mat_B = math.sqrt(occ_B)  # b|n> = √n|n-1>
        else:
            # Invalid transition
            return False

        matrix_element_product = mat_A * mat_B
        time_factor = self.hamiltonian.t * self.beta
        neighbors = self.lattice.get_neighbors(site)
        z = len(neighbors)

        # Reverse of insert_kink proposal ratio
        proposal_ratio = 1.0 / (2.0 * z)

        acceptance_ratio = proposal_ratio / (matrix_element_product * time_factor)
        acceptance_prob = min(1., acceptance_ratio)

        if acceptance_prob < self.rng.random():
            return False

        # Remove the hop events
        self.config.remove_element(site, hop_index)
        self.config.remove_element(neighbor_site, neighbor_hop_index)
        # Find and remove the worm event at current worm site
        worm_index = self._find_event_index(site, current_time, type_filter=worm_type)
        if worm_index is None:
            worm_index = self._find_event_index(site, current_time,
                                        eps=10*EPSILON, type_filter=worm_type)
            if worm_index is None:
                return False

        self.config.remove_element(site, worm_index)

        # ... remove hop events at site and neighbor_site ...

        # Reinsert the worm at the ORIGINAL site (before kink insertion)
        if worm_type == TYPE_WORM_HEAD:
            # Worm was at current site, move back to original site
            original_site = neighbor_site
            index_worm = self.config.insert_element(
                original_site, current_time,
                TYPE_WORM_HEAD, final_occ_B, occ_B,
                linked_site=original_site
            )
            self.config.worm_head_site = original_site
            self.config.worm_head_time = current_time
            self.config.worm_head_wpm = 0
        else:  # TYPE_WORM_TAIL
            original_site = neighbor_site
            index_worm = self.config.insert_element(
                original_site, current_time,
                TYPE_WORM_TAIL, final_occ_B, occ_B,
                linked_site=original_site
            )
            self.config.worm_tail_site = original_site
            self.config.worm_tail_time = current_time
            self.config.worm_tail_wpm = 0
        """
        # Find and remove the worm event at neighbor site
        worm_index = self._find_event_index(neighbor_site, current_time,
                                       type_filter=worm_type)
        if worm_index is None:
            worm_index = self._find_event_index(neighbor_site, current_time,
                                           eps=10*EPSILON, type_filter=worm_type)
            if worm_index is None:
                return False

        self.config.remove_element(neighbor_site, worm_index)

        # Reinsert the worm at the original site
        # Determine which site gets the worm back based on worm type
        if worm_type == TYPE_WORM_HEAD:
            # Worm was at neighbor site, move back to original site
            original_site = site
            # Occupations: we go from final_occ_B back to occ_B
            index_worm = self.config.insert_element(original_site, current_time,
                                               TYPE_WORM_HEAD, final_occ_B, occ_B,
                                               linked_site=original_site)
            self.config.worm_head_site = original_site
            self.config.worm_head_time = current_time
            self.config.worm_head_wpm = 0
        else:  # TYPE_WORM_TAIL
            original_site = site
            index_worm = self.config.insert_element(original_site, current_time,
                                               TYPE_WORM_TAIL, final_occ_B, occ_B,
                                               linked_site=original_site)
            self.config.worm_tail_site = original_site
            self.config.worm_tail_time = current_time
            self.config.worm_tail_wpm = 0
        """
        self.stats['deletekink_accepts'] += 1
        self._log(f"[DeleteKink] {'head' if move_head else 'tail'}; site {neighbor_site} --> {original_site}, time {current_time:.6f}")
        return True

    def monte_carlo_sweep(self, updates_per_sweep = 100):

        for _ in range(updates_per_sweep):
            if self.config.in_z_sector:
                self.insert_worm()
            else:
                u = self.rng.random()
                if u < self.move_prob:
                    self.move_worm()
                elif u < self.move_prob + self.glue_prob:
                    self.glue_worm()
                elif u < self.move_prob + self.glue_prob + self.kink_prob / 2:
                    self.insert_kink()
                else:
                    self.delete_kink()

        self.stats['sweeps'] += 1

    def get_acceptance_rates(self):
        rates = {}
        for key in ('insert', 'glue', 'move', 'insertkink', 'deletekink'):
            attempts = self.stats.get(f'{key}_attempts', 0)
            accepts = self.stats.get(f'{key}_accepts', 0)
            rates[key] = accepts / max(1, attempts)

        return rates

