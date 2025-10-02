# Example Usage Documentation

## New Fade Parameter Functionality

You can now set `fade_in_time` and `fade_in_and_out` parameters individually for each `Trace` instance. If not specified (set to `None`), the trace will use the default values from `make_plot()`.

### Examples:

```python
# Example 1: Trace with custom fade parameters
fast_fade_trace = Trace(
    data_df, 
    data_name='Fast Fade Cluster', 
    fade_in_time=2.0,  # Custom fade time of 2 Myr
    fade_in_and_out=True,  # Enable fade in and out
    min_size=5, 
    max_size=15, 
    color='red'
)

# Example 2: Trace with only custom fade_in_time
custom_fade_trace = Trace(
    data_df, 
    data_name='Custom Fade Cluster', 
    fade_in_time=10.0,  # Custom fade time of 10 Myr
    fade_in_and_out=None,  # Will use default from make_plot
    min_size=3, 
    max_size=12, 
    color='blue'
)

# Example 3: Trace using defaults (existing behavior)
default_trace = Trace(
    data_df, 
    data_name='Default Cluster', 
    # fade_in_time and fade_in_and_out not specified (None)
    min_size=2, 
    max_size=8, 
    color='green'
)

# Create collection with mixed fade settings
collection = TraceCollection([fast_fade_trace, custom_fade_trace, default_trace])

# When calling make_plot, the default values are used for traces with None values
plot_3d = Animate3D(data_collection=collection, figure_theme='dark')
fig = plot_3d.make_plot(
    time=np.arange(0, -50, -1),
    fade_in_time=5.0,  # Default for traces with fade_in_time=None
    fade_in_and_out=False,  # Default for traces with fade_in_and_out=None
    show=True
)
```

### Result:
- `fast_fade_trace`: Uses fade_in_time=2.0, fade_in_and_out=True
- `custom_fade_trace`: Uses fade_in_time=10.0, fade_in_and_out=False (from make_plot default)  
- `default_trace`: Uses fade_in_time=5.0, fade_in_and_out=False (both from make_plot defaults)