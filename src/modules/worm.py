import numpy as np
import bisect

from configuration import WormConfiguration
from utils import WormUtils

# Event types 

TYPE_HOP = 0
TYPE_WORM_HEAD = 1
TYPE_WORM_TAIL = 2
TYPE_DUMMY = 3

class WormAlgorithm:

    def __init__(
            self, 
            lattice,
            hamiltonian,
            beta,
            n_max = 50
            c_worm = 1.0
            e_off = 1.0
            seed = None
            epsilon_time = 1.0e-3
            glue_prob = 0.4
            move_prob = 0.6
            verbose = False
            ):

        self.lattece = lattice
        self.hamiltonian = hamiltonian
        self.beta = beta
        self.n_max = n_max
        self.c_worm = c_worm
        self.e_off = e_off

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
    
    self.stats['insert_attmepts'] += 1

    site = self.rng.integers(0, self.lattice.get_nsites())
    time = self.rng.uniform(0.0, self.beta)
    occ = self.config.get_occupation_at_time(site, time)
    new_occ = np.copy(occ)

    create_first = self.rng.random() < 0.5
    if create_first:
        new_occ += 1
        if new_occ > elf.n_max:
            return False
    else :
        if occ = 0:
            return False
        
        new_occ -= 1

    mat_prod = self._matrix_prod_from_occ_change(occ, new_occ)
    
    if mat_prod <= 0.0:
        return False

    acceptance_ratio = 2 * self.c_worm  * mat_prod
    acceptance_prob = min(1.0, acceptance_ratio)
    
    if acceptance_prob < self.rng.random():
        return False

    epsilon = max(1.0e-12 * self.beta, 1.0e-15 )
    t_tail = time
    t_head = time + epsilon
    
    index_tail = self.config.insert_element(site, t_tail, TYPE_WORM_TAIL, occ, new_occ, linked_site=site)
    index_head = self.config.insert_element(site, t_head, TYPE_WORM_HEAD, new_occ, occ, linked_site=site)

    self.config.worm_tail_site = site
    self.config.worm_tail_time = self.config.events[site][index_tail]['time']
    self.config.worm_head_site = site
    self.config.worm_head_time = self.config.events[site][index_head]['time']
    self.config.in_z_sector = False
    self.stats['insert_accepts'] += 1
    self._log(f'[Insert Worm] site = {site}; times = ({t_tail:.6f}, {t_head:.6f});  occ:{occ} --> {new_occ}')
    return True

    def glue_worm(self):
        """
        Cierra el worm si los extremo están abiertos o el ordenamiento del tiempo imaginario está invertido.
        """

        if self.config.in_z_sector:
            return False

        self.stats['glue_attemps'] += 1

        head_site = self.config.worm_head_site
        tail_site = self.config.worm_tail_site
        if head_site != tail_site:
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

        if not (close_enough or crossed):
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
        self.config.in_z_sector = True
        self.stats('glue_accepts') += 1
        self._log(f'[Glue] site = {site}; times = ({time_head:.6f},{time_tail:.6f}); occ: {occ} --> {occ_after}')
        return True


    def move_worm(self):

        eps = 1.0e-12
        if self.config.in_z_sector:
            return False

        self.stats['move_attempts'] += 1

        site = self.config.worm_head_site
        current_time = self.config.worm_head_time
        event_type = TYPE_WORM_HEAD if move_head else TYPE_WORM_TAIL

        prev_time = self.config.find_prev_event_time(site, current_time)
        next_time = self.config.find_next_event_time(site, current_time)

        E_before = self.config.compute_local_energy(site, current_time - eps)
        E_after = self.config.compute_local_energy(site, current_time + eps)
        delta_E = E_after - E_before

        rate = max( eps, self.e_off + np.abs(delta_E) )
        
        delta = WormUtils.exponential_deviate(rate, self.rng)
        proposed_time = current_time + delta if (self.rng.random() < 0.5) else current_time - delta
        if proposed_time <= prev_time or proposed_time >= next_time:
            return False

        E_new_before = self.config.compute_local_energy(site, proposed_time - eps)
        E_new_after = self.config.compute_local_energy(site, proposed_time + eps)

        E_new = 0.5 * (E_new_before + E_new_after)
        E_old = 0.5 * (E_before + E_after)
        delta_tau = proposed_time - current_time
        deltaS = (E_new - E_old)*delta_tau 

        acceptance_ratio = 0.0 if (deltaS > 700) else np.exp(-deltaS)
        acceptance_prob = min(1.0, acceptance_ratio)

        if (acceptance_prob < self.rng.random() ):
            return False

        index = self._find_event_index(site, current_time, type_filter = event_type)
        
        if index is None:
            return False
        event = self.config.events[site][index].copy()
        self.config.remove_element(site, index)
        new_index = self.config.insert_element(site, proposed_time, event['type'], event['occ_left'], event['occ_right'], event['linked_site'])

        if self.rng.random() < 0.5:
            self.config.worm_head_time = self.config.events[site][new_index]['time']
        else:
            self.config.worm_tail_time = self.config.events[site][new_index]['time']

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
 

        


