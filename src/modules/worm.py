import numpy as np
import bisect

from configuration import WormConfiguration
#from utils import WormUtils

# Event types 

TYPE_HOP = 0
TYPE_WORM_HEAD = 1
TYPE_WORM_TAIL = 2
TYPE_DUMMY = 3

EPSILON = 1.0e-15

class WormAlgorithm:

    def __init__(
            self, 
            lattice,
            hamiltonian,
            beta,
            n_max = 50
            c_worm = 1.0
            energy_off = 1.0
            seed = None
            epsilon_time = 1.0e-3
            glue_prob = 0.4
            move_prob = 0.6
            verbose = False
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
                'insert_attemps': 0,
                'insert_accepts': 0,
                'glue_attemps'  : 0,
                'glue_accepts'  : 0,
                'move_attempts' : 0,
                'move_accepts'  : 0,
                'sweeps'        : 0,
                }

        self.epsilon_time = epsilon_time * self.beta 

        self.glue_prob = glue_prob
        self.move_prob = move_prob

        self.verbose = verbose

        # helpers

        def _log(self, *args, **kwargs):
            if self.verbose:
                print(*args, **kwargs)

        def _find_event_index(self, site, time, type_filter = None):
            """
            Encontrar el indice siguiente que coincide con time dentro de la tolerancia
            """
            times = [e['time'] for e in self.config.events[site]]
            if not times:
                return None

            tnorm = self.config._norm_time(time)

            i = bisect.bisect_left(times, tnorm)

            for index in (i - 1, i, i + 1):
                if 0 <= index < len(times):
                    event = self.config.events[site]


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
        new_occ = np.copy(occ)                                 #We add change the occupation later 

        create_first = self.rng.random() < 0.5                 #We decide with equal probability whether we create of annihilate first
        if create_first:
            new_occ += 1
            if new_occ > self.n_max:
                return False                                   #Occupation out of boundaries (superior bound surpassed), we reject
        else :
            if occ = 0:                                        #Occupation out of boundaries (we can't annhilate particles if they do not exist), we reject 
                return False

            new_occ -= 1

        mat_prod = self._matrix_prod_from_occ_change(occ, new_occ) #|< occ | b_A | new_occ >|^2
    
        if mat_prod <= 0.0:                                    #It shouldn't be less equal zero, but just in case
            return False

        acceptance_ratio = 2 * self.c_worm  * mat_prod         #Acceptance ratio
        acceptance_prob = min(1.0, acceptance_ratio)           #Acceptance probability
    
        if acceptance_prob < self.rng.random():                #if the probability is less than a random number, we reject
            return False

        epsilon = max(1.0e-15 * self.beta, 1.0e-15 )           #We insert the worms a infinitesimal time appart
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
        self._log(f'[Insert Worm] site = {site}; times = ({t_tail:.6f}, {t_head:.6f});  occ:{occ} --> {new_occ}')
        return True

    def glue_worm(self):
        """
        Cierra el worm si los extremo están abiertos o el ordenamiento del tiempo imaginario está invertido.
        """

        if self.config.in_z_sector: #Reject if we are not in the G sector
            return False

        self.stats['glue_attemps'] += 1

        head_site = self.config.worm_head_site         
        tail_site = self.config.worm_tail_site
        if head_site != tail_site:            #tail and head must be in the same site
            return False

        site = head_site

        time_head = self.config._norm_time(self.config.worm_head_time)
        time_tail = self.config._nomr_time(self.config.worm_tail_time)

        diff = time_head - time_tail
        if diff > 0.5 * self.beta:
            diff -= self.beta
        elif diff <= 0.5 * self.beta
            diff += self.beta 

        close_enough = np.abs(diff) <= self.epsilon_time
        crossed = np.abs(diff) > (0.5 * self.beta - self.epsilon_time)

        if not (close_enough or crossed):         #Time must be the same within tolerance if not reject
            return False

        index_head = self._find_event_index(site, time_head, TYPE_WORM_HEAD)  
        index_tail = self._find_event_index(site, time_tail, TYPE_WORM_TAIL)
        if index_head is None or index_tail is None:
            return False

        ref = min(time_head, time_tail)

        occ = self.config.get_occupation_at_time(site, ref - 1.0e-12)
        new_occ = self.config.get_occupation_at_time(site, ref + 1.0e-12)

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
        self.stats('glue_accepts') += 1
        self._log(f'[Glue] site = {site}; times = ({time_head:.6f},{time_tail:.6f}); occ: {occ} --> {occ_after}')
        return True


    def move_worm(self, eps = EPSILON):
        
        tolerance = eps * self.beta
        if self.config.in_z_sector:           #Must be in the G sector
            return False

        self.stats['move_attempts'] += 1
        
        move_head = self.rng.random() < 0.5        #we decide with equal probability which worm we'll move
        forward = self.rng.rando() < 0.5
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
                pexp = -np.log( self.rng.random() ) /rate
                r_denom = rate if tau > pexp else 1.0
                r_num = E_R - E_L + self.energy_off if current_wpm == 0 else 1.0
            else:
                rate = E_L - E_R + self.energy_off
                pexp = -np.log( self.rng.random() ) / rate
                r_denom = rate if tau > pexp else 1.0
                r_num = self.energy_off if current_wpm ==  0 else 1.0
            proposed_time = current_time + min(pexp, tau)
        else: 
            tau = self.config.find_time_to_prev_element(site,current_time, include_neighbors = True)
            if E_R > E_L:
                rate = E_R + E_L - self.energy_off
                pexp = -np.log( self.rng.random() ) /rate
                r_denom = rate if tau > pexp else 1.0
                r_num = self.energy_off if current_wpm == 0 else 1.0
            else:
                rate = self.energy_off
                pexp = -np.log( self.rng.random() ) / rate
                r_denom = rate if tau > pexp else 1.0
                r_num = E_L - E_R + self.energy_off if current_wpm == 0 else 1.0
            proposed_time = -min(pexp, tau)

        proposed_time = self.config._norm_time(proposed_time)

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

    def monte_carlo_sweep(self, updates_per_sweep = 100):

        for _ in range(updates_per_swepp):
            if self.config.in_z_sector:
                self.insert_worm()
            else:
                u = self.rng.random()
                if u < self.move_prob:
                    self.move_worm()
                elif u < self.move_prob + self.glue_prob:
                    self.glue_worm()
                else:
                    self.insert_worm()

            self.stats['sweeps'] += 1

    def get_acceptance_rates(self):
        rates = {}
        for key in ('insert', 'glue', 'move'):
            attempts = self.stats.get(f'{key}_attempts', 0)
            accepts = self.stats.get(f'{key}_accepts', 0)
            rates[key] = accepts / max(1, attempts)

        return rates
 

        


