# ============================================================================
#
# This file is part of: AMAZED
#
# Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
#
# https://www.lam.fr/
#
# This software is a computer program whose purpose is to estimate the
# spectrocopic redshift of astronomical sources (galaxy/quasar/star)
# from there 1D spectrum.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
# ============================================================================
from contextlib import contextmanager

from pylibamazed.redshift import (CFlagWarning, CLog, CProcessFlowContext,
                                  ErrorCode, GlobalException, ScopeType)

zflag = CFlagWarning.GetInstance()
zlog = CLog.GetInstance()
ProcessFlowContext = CProcessFlowContext.GetInstance()


scopeStack = ProcessFlowContext.m_ScopeStack


def get_scope_level():
    return scopeStack.size()


def get_scope_type(scopetype):
    return scopeStack.get_type_value(scopetype.value)


def get_scope_SpectrumModel():
    try:
        spectrum_model = ProcessFlowContext.GetCurrentCategory()
    except GlobalException as e:
        errCode = ErrorCode(e.getErrorCode())
        if (errCode.value != ErrorCode.SCOPESTACK_ERROR):
            raise e from None
        spectrum_model = None

    # handle the classification stage, till empy string is not supported in C++ scope and resultstore
    if spectrum_model == "classification":
        spectrum_model = None
    return spectrum_model


def get_scope_stage():
    try:
        stage = ProcessFlowContext.GetCurrentStage()
    except GlobalException as e:
        errCode = ErrorCode(e.getErrorCode())
        if (errCode.value != ErrorCode.SCOPESTACK_ERROR):
            raise e from None
        stage = None

    # handle the init stage (root scope)
    if get_scope_SpectrumModel() is None and stage is None:
        stage = "init"
    return stage


def log_scope(ending=False):
    scope_type = scopeStack.get_current_type()
    if scope_type == ScopeType.UNDEFINED:
        return
    scope = scopeStack.back()
    prefix = "Starting" if not ending else "Ending"
    message = f"{prefix} {scope_type.name}: {scope}"
    zlog.LogInfo(message)


@contextmanager
def push_scope(scope, scope_type=None):
    if scope_type is None:
        scopeStack.push_back(scope)
    else:
        scopeStack.push_back(scope, scope_type.value)
        log_scope()
    try:
        yield
    finally:
        if scope_type is not None:
            log_scope(ending=True)
        scopeStack.pop_back()
