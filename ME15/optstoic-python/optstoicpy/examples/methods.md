# optStoic Python tutorial
- Author: Chiam Yu Ng
- Reference:  Chiam Yu Ng, Lin Wang, Anupam Chowdhury and Costas D. Maranas, 2018 (under review)
- Details on the algorithms is provided in the paper.

## Methods
1. Setup your optStoic `design equation`.

2. Create a `user_defined_export_rxns_Sji` dictionary of exchange reactions for each metabolite that participate in the equation. E.g., 
```python
user_defined_export_rxns_Sji = {
    'EX_glc': {'C00031': -1.0},
    'EX_nad': {'C00003': -1.0},
    'EX_adp': {'C00008': -1.0},
    'EX_phosphate': {'C00009': -1.0},
    'EX_pyruvate': {'C00022': -1.0},
    'EX_nadh': {'C00004': -1.0},
    'EX_atp': {'C00002': -1.0},
    'EX_h2o': {'C00001': -1.0},
    'EX_hplus': {'C00080': -1.0},
    'EX_nadp': {'C00006': -1.0},
    'EX_nadph': {'C00005': -1.0}
    }
```

3. Create a `specific_bounds` dictionary for the upper and lower bound to the exchange flux. 
```python
specific_bounds = {'EX_glc': {'LB': -1, 'UB': -1},
                'EX_pyruvate': {'LB': 2, 'UB': 2},
                'EX_nad': {'LB': -2, 'UB': 0},
                'EX_nadh': {'LB': 0, 'UB': 2},
                'EX_nadp': {'LB': -2, 'UB': 0},
                'EX_nadph': {'LB': 0, 'UB': 2},
                'EX_adp': {'LB': -1, 'UB': -1},
                'EX_phosphate': {'LB': -1, 'UB': -1},
                'EX_atp': {'LB': 1, 'UB': 1},
                'EX_h2o': {'LB': 1, 'UB': 1},
                'EX_hplus': {'LB': -10, 'UB': 10}}
```

4. Setup your custom constraints for the pathway.
```
v('EX_nadph') + v('EX_nadh') = 2;
v('EX_nadp') + v('EX_nad') = -2;
v('EX_nadh') + v('EX_nad') = 0;
v('EX_nadph') + v('EX_nadp') = 0;
```
can be written as:
```python    
custom_flux_constraints = [
        {'constraint_name': 'nadphcons1',
         'reactions': ['EX_nadph', 'EX_nadh'],
         'UB': 2,
         'LB': 2},
        {'constraint_name': 'nadphcons2',
        'reactions': ['EX_nadp', 'EX_nad'],
        'UB': -2,
        'LB': -2},
        {'constraint_name': 'nadphcons3',
        'reactions': ['EX_nadh', 'EX_nad'],
        'UB': 0,
        'LB': 0},
        {'constraint_name': 'nadphcons4',
        'reactions': ['EX_nadph', 'EX_nadp'],
        'UB': 0,
        'LB': 0}]
```

5. Perform blocked reaction analysis (see `optstoicpy.script.data_preprocessing.blocked_reactions_analysis` details).

6. Perform internal loop analysis to obtain the `Nint` matrix (see `optstoicpy.script.data_preprocessing.test_internal_loop_analysis` for details).

7. Setup optStoic database loading (`optstoicpy.core.database.load_db_v3`).

8. Finally, setup the optStoic simulation (`optstoicpy.script.opstoic`).

9. Visualize pathways (`optstoicpy.core.drawpathway`) and analyze pathways (`optstoicpy.script.pathway_analysis`)
