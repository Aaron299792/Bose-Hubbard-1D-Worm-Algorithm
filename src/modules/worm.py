import bisect
import math
import random

from configuration import WormConfiguration
from configuration import TYPE_HOP, TYPE_WORM_TAIL, TYPE_WORM_HEAD, TYPE_WORM_DUMMY, EPSILON

class WormAlgorithm:

    def __init__(
            self,
            lattice,
            hamiltonian,
            beta,
            n_max = 50,
            c_worm = 1.0,
            energy_off = 1.0,
            seed = None,
            epsilon_time = 1.0e-3,
            kink_prob = 0.25,
            glue_prob = 0.4,
            move_prob = 0.35,
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
        """
        Cierra el worm si los extremo están abiertos o el ordenamiento del tiempo imaginario está invertido.
        """

        if self.config.in_z_sector: #Reject if we are not in the G sector
            return False

        self.stats['glue_attempts'] += 1

        head_site = self.config.worm_head_site
        tail_site = self.config.worm_tail_site
        if head_site != tail_site:            #tail and head must be in the same site
            return False

        site = head_site

        time_head = self.config._norm_time(self.config.worm_head_time)
        time_tail = self.config._norm_time(self.config.worm_tail_time)

        time_diff = self.config._time_distance(time_tail, time_head)
        if time_diff > self.epsilon_time:
            return False

        index_head = self._find_event_index(site, time_head, type_filter = TYPE_WORM_HEAD)
        index_tail = self._find_event_index(site, time_tail, type_filter = TYPE_WORM_TAIL)
        if index_head is None or index_tail is None:
            return False

        ref = min(time_head, time_tail)

        occ = self.config.get_occupation_at_time(site, ref - 1000*eps)
        new_occ = self.config.get_occupation_at_time(site, ref + 1000*eps)

        mat_prod = self._matrix_prod_from_occ_change(occ, new_occ)
        if mat_prod <= 0.0:
            return False

        acceptance_ratio = 1.0 / (2.0 * self.c_worm * mat_prod)
        acceptance_prob = min(1.0, acceptance_ratio)

        if acceptance_prob < self.rng.random():
            return False

        for index in sorted((index_head, index_tail), reverse=True):
            self.config.remove_element(site, index)

        self.config.worm_head_site = -1
        self.config.worm_head_time = -1.0
        self.config.worm_tail_site = -1
        self.config.worm_tail_time = -1.0
        self.config.worm_head_wpm = 0
        self.config.worm_tail_wpm = 0
        self.config.in_z_sector = True
        self.stats['glue_accepts'] += 1
        self._log(f'[Glue] site = {site}; times = ({time_head:.6f},{time_tail:.6f}); occ: {occ} --> {new_occ}')
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
        site, current_time, current_wpm, worm_type = 0, 0.0, 0, 0

        move_head = self.rng.random() < 0.5

        if move_head:
            site =self.config.worm_head_site
            current_time = self.config.worm_head_time
            current_wpm  = self.config.worm_head_wpm
            worm_type = TYPE_WORM_HEAD
        else:
            site = self.config.worm_tail_site
            current_time = self.config.worm_tail_time
            current_wpm  = self.config.worm_tail_wpm
            worm_type = TYPE_WORM_TAIL

        if current_wpm != 0:
            return False

        neighbors = self.lattice.get_neighbors(site)
        if not neighbors:
            return False

        neighbor_site = self.rng.choice(neighbors)
        occ_A = self.config.get_occupation_at_time(site, current_time)
        occ_B = self.config.get_occupation_at_time(neighbor_site, current_time)

        final_occ_A, final_occ_B, mat_A, mat_B = 0, 0, 0, 0

        if worm_type == TYPE_WORM_HEAD:
            final_occ_A = occ_A - 1 #particle number conservation
            final_occ_B = occ_B + 1
            mat_A = math.sqrt(occ_A)
            mat_B = math.sqrt(occ_B + 1)

        else:
            final_occ_A = occ_A + 1
            final_occ_B = occ_B - 1
            mat_A = math.sqrt(occ_A + 1)
            mat_B = math.sqrt(occ_B)

        if (final_occ_A < 0 or final_occ_A > self.n_max or
            final_occ_B < 0 or final_occ_B > self.n_max):
            return False

        matrix_element_product = mat_A * mat_B
        time_factor = self.hamiltonian.t* self.beta
        z = len(neighbors)
        proposal_ratio = 2.0*z

        acceptance_ratio = matrix_element_product * time_factor * proposal_ratio
        acceptance_prob = min(1., acceptance_ratio)

        if acceptance_prob < self.rng.random():
            return False

        index_old = self._find_event_index(site, current_time, type_filter=worm_type)
        if index_old is None:
            return False

        old_event = self.config.events[site][index_old].copy()
        self.config.remove_element(site, index_old)

        index_hop_A = self.config.insert_element(site, current_time, TYPE_HOP, occ_A, final_occ_A, linked_site= neighbor_site)
        index_hop_B = self.config.insert_element(neighbor_site, current_time, TYPE_HOP, occ_B, final_occ_B, linked_site = site)

        epsilon = max(EPSILON * self.beta, EPSILON)

        place_after = self.rng.random() < 0.5
        if place_after:
            worm_time = current_time + epsilon
            new_wpm = 1
        else:
            worm_time = current_time - epsilon
            new_wpm = -1

        index_worm = self.config.insert_element(neighbor_site, worm_time, worm_type, final_occ_B, occ_B, linked_site=neighbor_site)

        if move_head:
            self.config.worm_head_site = neighbor_site
            self.config.worm_head_time = self.config.events[neighbor_site][index_worm]['time']
            self.config.worm_head_wpm = new_wpm
        else:
            self.config.worm_tail_site = neighbor_site
            self.config.worm_tail_time = self.config.events[neighbor_site][index_worm]['time']

            self.config.worm_tail_wpm = new_wpm


        self.stats['insertkink_accepts'] += 1
        self._log(f"[InsertKink] {'head' if move_head else 'tail'}; site {site} --> {neighbor_site}, time {current_time:.6f}")
        return True

    def delete_kink(self):
        if self.config.in_z_sector:
            return False

        self.stats['deletekink_attempts'] += 1
        site, current_time, current_wpm, worm_type = 0, 0.0, 0, 0

        move_head = self.rng.random() < 0.5

        if move_head:
            site =self.config.worm_head_site
            current_time = self.config.worm_head_time
            current_wpm  = self.config.worm_head_wpm
            worm_type = TYPE_WORM_HEAD
        else:
            site = self.config.worm_tail_site
            current_time = self.config.worm_tail_time
            current_wpm  = self.config.worm_tail_wpm
            worm_type = TYPE_WORM_TAIL

        if current_wpm != 0:
            return False

        hop_index = self._find_event_index(site, current_time, type_filter=TYPE_HOP)
        if hop_index is None:
            return False

        hop_event = self.config.events[site][hop_index]
        neighbor_site = hop_event['linked_site']

        neighbor_hop_index = self._find_event_index(neighbor_site, current_time, type_filter=TYPE_HOP)
        if neighbor_hop_index is None:
            return False

        neighbor_hop_event = self.config.events[neighbor_site][neighbor_hop_index]

        if neighbor_hop_event['linked_site'] != site:
            return False

        occ_A = hop_event['occ_left']
        final_occ_A = hop_event['occ_right']
        occ_B = neighbor_hop_event['occ_left']
        final_occ_B = neighbor_hop_event['occ_right']

        if worm_type == TYPE_WORM_HEAD:
            mat_A = math.sqrt(occ_A)
            mat_B = math.sqrt(occ_B + 1)
        else:
            mat_A = math.sqrt(occ_A + 1)
            mat_B = math.sqrt(occ_B)

        matrix_element_product = mat_A * mat_B
        time_factor = self.hamiltonian.t* self.beta
        neighbors = self.lattice.get_neighbors(site)
        z = len(neighbors)
        proposal_ratio = 2.0 / float(z)

        acceptance_ratio = proposal_ratio / (matrix_element_product * time_factor)
        acceptance_prob = min(1., acceptance_ratio)

        if acceptance_prob < self.rng.random():
            return False

        self.config.remove_element(site, hop_index)
        self.config.remove_element(neighbor_site, neighbor_hop_index)

        worm_index = self._find_event_index(neighbor_site, current_time, type_filter=worm_type)
        if worm_index is None:
            return False

        self.config.remove_element(neighbor_site, worm_index)

        original_site = neighbor_site if move_head else site
        new_worm_time = current_time

        if worm_type == TYPE_WORM_HEAD:
            index_worm = self.config.insert_element(original_site, new_worm_time,
                                                    TYPE_WORM_HEAD, final_occ_B, occ_B,
                                                    linked_site= original_site)
            self.config.worm_head_site = original_site
            self.config.worm_head_time = new_worm_time
            self.config.worm_head_wpm = 0

        else:
            index_worm = self.config.insert_element(original_site, new_worm_time,
                                                    TYPE_WORM_TAIL, final_occ_B, occ_B,
                                                    linked_site=original_site)
            self.config.worm_tail_site = original_site
            self.config.worm_tail_time = new_worm_time
            self.config.worm_tail_wpm = 0

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
                elif u < self.move_prob + self.glue_prob + self.kink_prob:
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
