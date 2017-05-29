from time import time
import logging
import numpy as np
import os
import abc
import sys


class GenericSampler(object):
    __metaclass__ = abc.ABCMeta

    def fit(self, kwargs):
        """ Fits a given model using the Sampler.

        Parameters
        ----------
        kwargs : dict
            Dictionary of keyword arguments utilised by the fitters

        Returns
        -------
        dict
            A dictionary of results containing:
                - *chain*: the chain
                - *weights*: chain weights if applicable
        """
        raise NotImplementedError()


class EnsembleSampler(GenericSampler):
    def __init__(self, num_walkers=None, num_steps=5000, num_burn=2000,
                 temp_dir=None, save_interval=300):
        """ Uses ``emcee`` and the `EnsembleSampler
        <http://dan.iel.fm/emcee/current/api/#emcee.EnsembleSampler>`_ to fit the supplied
        model.

        This method sets an emcee run using the ``EnsembleSampler`` and manual
        chain management to allow for low to medium dimensional models. MPI running
        is detected automatically for less hassle, and chain progress is serialised
        to disk automatically for convenience.

        Parameters
        ----------
        num_walkers : int, optional
            The number of walkers to run. If not supplied, it defaults to eight times the
            framework dimensionality
        num_steps : int, optional
            The number of steps to run
        num_burn : int, optional
            The number of steps to discard for burn in
        temp_dir : str
            If set, specifies a directory in which to save temporary results, like the emcee chain
        save_interval : float
            The amount of seconds between saving the chain to file. Setting to ``None``
            disables serialisation.
        """

        self.logger = logging.getLogger(__name__)
        import emcee
        self.chain = None
        self.pool = None
        self.master = True
        self.num_steps = num_steps
        self.num_burn = num_burn
        self.temp_dir = temp_dir
        if temp_dir is not None and not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        self.save_interval = save_interval
        self.num_walkers = num_walkers

    def fit(self, kwargs):
        """ Runs the sampler over the model and returns the flat chain of results

        Parameters
        ----------
        kwargs : dict
            Containing the following information at a minimum:

            - log_posterior : function
                A function which takes a list of parameters and returns
                the log posterior
            - start : function|list|ndarray
                Either a starting position, or a function that can be called
                to generate a starting position
            - save_dims : int, optional
                Only return values for the first ``save_dims`` parameters.
                Useful to remove numerous marginalisation parameters if running
                low on memory or hard drive space.
            - uid : str, optional
                A unique identifier used to differentiate different fits
                if two fits both serialise their chains and use the
                same temporary directory
        Returns
        -------
        dict
            A dictionary with key "chains" containing the final
            flattened chain of dimensions
             ``(num_dimensions, num_walkers * (num_steps - num_burn))``
        """
        log_posterior = kwargs.get("log_posterior")
        start = kwargs.get("start")
        save_dims = kwargs.get("save_dims")
        uid = kwargs.get("uid")
        assert log_posterior is not None
        assert start is not None
        from emcee.utils import MPIPool
        import emcee
        try:  # pragma: no cover
            self.pool = MPIPool()
            if not self.pool.is_master():
                self.logger.info("Slave waiting")
                self.master = False
                self.pool.wait()
                sys.exit(0)
            else:
                self.logger.info("MPIPool successful initialised and master found. "
                                 "Running with %d cores." % self.pool.size)
        except ImportError:
            self.logger.info("mpi4py is not installed or not configured properly. "
                             "Ignore if running through python, not mpirun")
        except ValueError as e:  # pragma: no cover
            self.logger.info("Unable to start MPI pool, expected normal python execution")
            self.logger.info(str(e))

        if callable(start):
            num_dim = start().size
        else:
            num_dim = start.size
        if self.num_walkers is None:
            self.num_walkers = num_dim * 4
            self.num_walkers = max(self.num_walkers, 20)

        self.logger.debug("Fitting framework with %d dimensions" % num_dim)

        self.logger.info("Using Ensemble Sampler")
        sampler = emcee.EnsembleSampler(self.num_walkers, num_dim,
                                        log_posterior,
                                        pool=self.pool, live_dangerously=True)

        emcee_wrapper = EmceeWrapper(sampler)
        flat_chain = emcee_wrapper.run_chain(self.num_steps, self.num_burn,
                                             self.num_walkers, num_dim,
                                             start=start,
                                             save_dim=save_dims,
                                             temp_dir=self.temp_dir,
                                             uid=uid,
                                             save_interval=self.save_interval)
        self.logger.debug("Fit finished")
        if self.pool is not None:  # pragma: no cover
            self.pool.close()
            self.logger.debug("Pool closed")

        return {"chain": flat_chain}


class EmceeWrapper(object):
    def __init__(self, sampler):
        self.sampler = sampler
        self.logger = logging.getLogger(__name__)
        self.chain = None

    def run_chain(self, num_steps, num_burn, num_walkers, num_dim, start=None, save_interval=300,
                  save_dim=None, temp_dir=None, uid="ensemble"):
        assert num_steps > num_burn, "num_steps has to be larger than num_burn"
        if save_dim is not None:
            assert save_dim <= num_dim, "You cannot save more dimensions than you actually have"
        else:
            save_dim = num_dim

        past_chain = None
        pos = None
        if temp_dir is not None:
            self.logger.debug("Looking in temp dir %s" % temp_dir)
            chain_file = temp_dir + os.sep + uid + "_ens.chain.npy"
            position_file = temp_dir + os.sep + uid + "_ens.pos.npy"
            try:
                pos = np.load(position_file)
                past_chain = np.load(chain_file)
                self.logger.info("Found chain of %d steps" % past_chain.shape[1])
            except IOError:
                self.logger.info("Prior chain and/or does not exist. Looked in %s" % position_file)

        if start is None and pos is None:
            raise ValueError("You need to have either a starting function or existing chains")

        if pos is None:
            pos = start(num_walkers=num_walkers)

        step = 0
        self.chain = np.zeros((num_walkers, num_steps, save_dim))
        if past_chain is not None:
            step = min(past_chain.shape[1], num_steps)
            num = num_steps - step
            self.chain[:, :step, :] = past_chain[:, :step, :]
            self.logger.debug("A further %d steps are required" % num)
        else:
            num = num_steps
            self.logger.debug("Running full chain of %d steps" % num)

        t = time()

        if step == num_steps:
            self.logger.debug("Returning serialised data from %s" % temp_dir)
            return self.get_results(num_burn)
        else:
            self.logger.debug("Starting sampling. Saving to %s ever %d seconds"
                              % (temp_dir, save_interval))

        for result in self.sampler.sample(pos, iterations=num, storechain=False):
            self.chain[:, step, :] = result[0][:, :save_dim]
            step += 1
            if step == 1 or temp_dir is not None and save_interval is not None:
                t2 = time()
                if temp_dir is not None and \
                        (step == 1 or t2 - t > save_interval or step == num_steps):
                    t = t2
                    position = result[0]
                    np.save(position_file, position)
                    np.save(chain_file, self.chain[:, :step, :])
                    self.logger.debug("Saving chain with %d steps" % step)
        return self.get_results(num_burn)

    def get_results(self, num_burn):
        return self.chain[:, num_burn:, :].reshape((-1, self.chain.shape[2]))
