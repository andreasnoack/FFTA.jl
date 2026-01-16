#!/usr/bin/env julia

"""
Generate HTML report with embedded Plotly.js charts for FFTA vs FFTW benchmarks
Uses client-side JavaScript to render interactive plots from JSON data
"""

# Ensure packages are installed
import Pkg
Pkg.instantiate()

using JSON
using Dates
using Primes

function generate_html_report()
    # Load 1D results
    ffta_file = joinpath(@__DIR__, "results_ffta.json")
    fftw_file = joinpath(@__DIR__, "results_fftw.json")

    if !isfile(ffta_file) || !isfile(fftw_file)
        error("1D Benchmark results not found. Run benchmarks first!")
    end

    ffta_results = JSON.parsefile(ffta_file)
    fftw_results = JSON.parsefile(fftw_file)

    # Load 2D results (optional)
    ffta_2d_file = joinpath(@__DIR__, "results_ffta_2d.json")
    fftw_2d_file = joinpath(@__DIR__, "results_fftw_2d.json")

    has_2d_results = isfile(ffta_2d_file) && isfile(fftw_2d_file)
    ffta_2d_results = has_2d_results ? JSON.parsefile(ffta_2d_file) : nothing
    fftw_2d_results = has_2d_results ? JSON.parsefile(fftw_2d_file) : nothing

    # Embed JSON data in JavaScript
    ffta_json = JSON.json(ffta_results)
    fftw_json = JSON.json(fftw_results)
    ffta_2d_json = has_2d_results ? JSON.json(ffta_2d_results) : "null"
    fftw_2d_json = has_2d_results ? JSON.json(fftw_2d_results) : "null"

    # Start HTML document
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>FFTA.jl vs FFTW.jl Performance Benchmark</title>
        <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
        <style>
            body {
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
                max-width: 1400px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f5f5f5;
                color: #333;
            }
            h1 {
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }
            h2 {
                color: #34495e;
                margin-top: 30px;
                border-bottom: 2px solid #95a5a6;
                padding-bottom: 5px;
            }
            h3 {
                color: #7f8c8d;
            }
            .summary {
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 20px 0;
            }
            .plot-container {
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 20px 0;
            }
            .plot {
                width: 100%;
                height: 600px;
            }
            .category-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(600px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                background-color: white;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            th, td {
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }
            th {
                background-color: #3498db;
                color: white;
                font-weight: 600;
            }
            tr:hover {
                background-color: #f5f5f5;
            }
            .faster {
                color: #27ae60;
                font-weight: 600;
            }
            .slower {
                color: #e74c3c;
                font-weight: 600;
            }
            .metadata {
                background-color: #ecf0f1;
                padding: 15px;
                border-radius: 4px;
                margin: 10px 0;
                font-family: monospace;
                font-size: 0.9em;
            }
            .info-box {
                background-color: #e8f4f8;
                border-left: 4px solid #3498db;
                padding: 15px;
                margin: 20px 0;
                border-radius: 4px;
            }
        </style>
    </head>
    <body>
        <h1>FFTA.jl vs FFTW.jl Performance Benchmark Report</h1>

        <div class="info-box">
            <strong>Note:</strong> This benchmark compares FFTA.jl (a pure Julia FFT implementation) against
            FFTW.jl (Julia bindings to the FFTW C library). Each package was benchmarked in a separate
            Julia process to ensure isolation and prevent interference.
        </div>

        <div class="summary">
            <h2>Benchmark Configuration</h2>
            <div class="metadata">
                <strong>Benchmarking Tool:</strong> BenchmarkTools.jl<br>
                <strong>Samples per size:</strong> 100<br>
                <strong>Evaluations per sample:</strong> 10<br>
                <strong>Data type:</strong> ComplexF64<br>
                <strong>Total array sizes tested:</strong> $(length(ffta_results["data"]))<br>
                <strong>Report generated:</strong> $(now())<br>
            </div>

            <h3>Array Size Categories</h3>
            <ul>
                <li><strong>Odd Powers of 2:</strong> 2¹, 2³, 2⁵, 2⁷, 2⁹, 2¹¹, 2¹³, 2¹⁵ (2, 8, 32, 128, 512, 2048, 8192, 32768)</li>
                <li><strong>Even Powers of 2:</strong> 2², 2⁴, 2⁶, 2⁸, 2¹⁰, 2¹², 2¹⁴ (4, 16, 64, 256, 1024, 4096, 16384)</li>
                <li><strong>Powers of 3:</strong> 3¹, 3², 3³, 3⁴, 3⁵, 3⁶, 3⁷, 3⁸, 3⁹ (3, 9, 27, 81, 243, 729, 2187, 6561, 19683)</li>
                <li><strong>Composite:</strong> 3, 12, 60, 300, 2100, 23100 (cumulative products of 3,4,5,5,7,11)</li>
                <li><strong>Primes:</strong> 20 logarithmically-spaced prime numbers up to 20,000</li>
            </ul>
        </div>

        <h2>1D FFT Performance Visualizations</h2>

        <div class="plot-container">
            <h3>Runtime/N vs N (All Categories)</h3>
            <div id="plot-combined" class="plot"></div>
        </div>

        <div class="plot-container">
            <h3>Absolute Runtime (All Categories)</h3>
            <div id="plot-absolute" class="plot"></div>
        </div>

        <div class="category-grid">
            <div class="plot-container">
                <h3>Odd Powers of 2</h3>
                <div id="plot-odd" class="plot" style="height: 400px;"></div>
            </div>
            <div class="plot-container">
                <h3>Even Powers of 2</h3>
                <div id="plot-even" class="plot" style="height: 400px;"></div>
            </div>
            <div class="plot-container">
                <h3>Powers of 3</h3>
                <div id="plot-power3" class="plot" style="height: 400px;"></div>
            </div>
            <div class="plot-container">
                <h3>Composite</h3>
                <div id="plot-composite" class="plot" style="height: 400px;"></div>
            </div>
            <div class="plot-container">
                <h3>Prime Numbers</h3>
                <div id="plot-primes" class="plot" style="height: 400px;"></div>
            </div>
        </div>

        <h2>1D FFT Detailed Results</h2>
        <div id="results-tables"></div>

        <div id="2d-section"></div>

        <div class="summary" style="margin-top: 40px;">
            <h2>Interpretation</h2>
            <p>
                The plots show <strong>Runtime/N vs N</strong>, which normalizes the runtime by the array size.
                For an optimal FFT implementation with O(N log N) complexity, this should scale as O(log N).
            </p>
            <p>
                <strong>Key Observations:</strong>
            </p>
            <ul>
                <li>Powers of 2 typically show the best performance for FFT algorithms due to the Cooley-Tukey radix-2 algorithm</li>
                <li>Powers of 3 may use different factorization strategies</li>
                <li>Composite numbers test the FFT implementation's ability to handle general factorizations</li>
                <li>FFTA.jl is a pure Julia implementation, while FFTW.jl wraps optimized C code</li>
                <li>Interactive plots allow zooming and hovering for detailed inspection</li>
            </ul>
        </div>

        <script>
            // Embedded benchmark data
            const fftaResults = $ffta_json;
            const fftwResults = $fftw_json;
            const ffta2dResults = $ffta_2d_json;
            const fftw2dResults = $fftw_2d_json;

            const categories = {
                'odd_power_of_2': { name: 'Odd Powers of 2', color: 'blue' },
                'even_power_of_2': { name: 'Even Powers of 2', color: 'red' },
                'power_of_3': { name: 'Powers of 3', color: 'green' },
                'composite': { name: 'Composite', color: 'purple' },
                'prime': { name: 'Prime Numbers', color: 'orange' }
            };

            const categories2d = {
                'square_power_of_2': { name: 'Square Powers of 2', color: 'blue' },
                'square_power_of_3': { name: 'Square Powers of 3', color: 'green' },
                'rect_power_of_2': { name: 'Rectangular Powers of 2', color: 'red' },
                'mixed_power2_power3': { name: 'Power of 2 × Power of 3', color: 'purple' },
                'mixed_power2_prime': { name: 'Power of 2 × Prime', color: 'orange' },
                'prime_prime': { name: 'Prime × Prime', color: 'brown' },
                'composite_2d': { name: 'Composite 2D', color: 'pink' }
            };

            // Helper function to filter data by category
            function filterByCategory(data, category) {
                return data.filter(d => d.category === category);
            }

            // Create combined Runtime/N plot
            function createCombinedPlot() {
                const traces = [];

                for (const [catKey, catInfo] of Object.entries(categories)) {
                    const fftaCat = filterByCategory(fftaResults.data, catKey);
                    const fftwCat = filterByCategory(fftwResults.data, catKey);

                    if (fftaCat.length > 0) {
                        traces.push({
                            x: fftaCat.map(d => d.size),
                            y: fftaCat.map(d => d.runtime_per_element * 1e9),
                            name: 'FFTA: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, color: catInfo.color },
                            line: { width: 2, color: catInfo.color }
                        });
                    }

                    if (fftwCat.length > 0) {
                        traces.push({
                            x: fftwCat.map(d => d.size),
                            y: fftwCat.map(d => d.runtime_per_element * 1e9),
                            name: 'FFTW: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, symbol: 'square', color: catInfo.color },
                            line: { width: 2, dash: 'dash', color: catInfo.color }
                        });
                    }
                }

                const layout = {
                    xaxis: { title: 'Array Length (N)', type: 'log' },
                    yaxis: { title: 'Runtime / N (nanoseconds)', type: 'log' },
                    hovermode: 'closest',
                    showlegend: true,
                    legend: { x: 1.05, y: 1 }
                };

                Plotly.newPlot('plot-combined', traces, layout, { responsive: true });
            }

            // Create absolute runtime plot
            function createAbsolutePlot() {
                const traces = [];

                for (const [catKey, catInfo] of Object.entries(categories)) {
                    const fftaCat = filterByCategory(fftaResults.data, catKey);
                    const fftwCat = filterByCategory(fftwResults.data, catKey);

                    if (fftaCat.length > 0) {
                        traces.push({
                            x: fftaCat.map(d => d.size),
                            y: fftaCat.map(d => d.median_time * 1e6),
                            name: 'FFTA: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, color: catInfo.color },
                            line: { width: 2, color: catInfo.color }
                        });
                    }

                    if (fftwCat.length > 0) {
                        traces.push({
                            x: fftwCat.map(d => d.size),
                            y: fftwCat.map(d => d.median_time * 1e6),
                            name: 'FFTW: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, symbol: 'square', color: catInfo.color },
                            line: { width: 2, dash: 'dash', color: catInfo.color }
                        });
                    }
                }

                const layout = {
                    xaxis: { title: 'Array Length (N)', type: 'log' },
                    yaxis: { title: 'Runtime (microseconds)', type: 'log' },
                    hovermode: 'closest',
                    showlegend: true,
                    legend: { x: 1.05, y: 1 }
                };

                Plotly.newPlot('plot-absolute', traces, layout, { responsive: true });
            }

            // Create category-specific plots
            function createCategoryPlot(catKey, divId) {
                const fftaCat = filterByCategory(fftaResults.data, catKey);
                const fftwCat = filterByCategory(fftwResults.data, catKey);

                if (fftaCat.length === 0 && fftwCat.length === 0) return;

                const traces = [];
                const catInfo = categories[catKey];

                if (fftaCat.length > 0) {
                    traces.push({
                        x: fftaCat.map(d => d.size),
                        y: fftaCat.map(d => d.runtime_per_element * 1e9),
                        name: 'FFTA.jl',
                        type: 'scatter',
                        mode: 'lines+markers',
                        marker: { size: 10, color: catInfo.color },
                        line: { width: 2, color: catInfo.color }
                    });
                }

                if (fftwCat.length > 0) {
                    traces.push({
                        x: fftwCat.map(d => d.size),
                        y: fftwCat.map(d => d.runtime_per_element * 1e9),
                        name: 'FFTW.jl',
                        type: 'scatter',
                        mode: 'lines+markers',
                        marker: { size: 10, symbol: 'square', color: 'gray' },
                        line: { width: 2, color: 'gray' }
                    });
                }

                const layout = {
                    xaxis: { title: 'Array Length (N)', type: 'log' },
                    yaxis: { title: 'Runtime / N (nanoseconds)', type: 'log' },
                    hovermode: 'closest',
                    showlegend: true
                };

                Plotly.newPlot(divId, traces, layout, { responsive: true });
            }

            // Create detailed results tables
            function createResultsTables() {
                let html = '';

                for (const [catKey, catInfo] of Object.entries(categories)) {
                    const fftaCat = filterByCategory(fftaResults.data, catKey);
                    const fftwCat = filterByCategory(fftwResults.data, catKey);

                    if (fftaCat.length === 0 && fftwCat.length === 0) continue;

                    html += '<h3>' + catInfo.name + '</h3>';
                    html += '<table>';
                    html += '<thead><tr>';
                    html += '<th>Array Size (N)</th>';
                    html += '<th>FFTA Time (μs)</th>';
                    html += '<th>FFTW Time (μs)</th>';
                    html += '<th>FFTA Runtime/N (ns)</th>';
                    html += '<th>FFTW Runtime/N (ns)</th>';
                    html += '<th>Speedup</th>';
                    html += '</tr></thead>';
                    html += '<tbody>';

                    const allSizes = [...new Set([...fftaCat.map(d => d.size), ...fftwCat.map(d => d.size)])].sort((a, b) => a - b);

                    for (const size of allSizes) {
                        const fftaData = fftaCat.find(d => d.size === size);
                        const fftwData = fftwCat.find(d => d.size === size);

                        if (fftaData && fftwData) {
                            const fftaTime = (fftaData.median_time * 1e6).toFixed(3);
                            const fftwTime = (fftwData.median_time * 1e6).toFixed(3);
                            const fftaPerN = (fftaData.runtime_per_element * 1e9).toFixed(3);
                            const fftwPerN = (fftwData.runtime_per_element * 1e9).toFixed(3);
                            const speedup = fftwData.median_time / fftaData.median_time;
                            const speedupClass = speedup > 1 ? 'faster' : 'slower';
                            const speedupText = speedup > 1
                                ? speedup.toFixed(2) + 'x (FFTA faster)'
                                : (1/speedup).toFixed(2) + 'x (FFTW faster)';

                            html += '<tr>';
                            html += '<td>' + size + '</td>';
                            html += '<td>' + fftaTime + '</td>';
                            html += '<td>' + fftwTime + '</td>';
                            html += '<td>' + fftaPerN + '</td>';
                            html += '<td>' + fftwPerN + '</td>';
                            html += '<td class="' + speedupClass + '">' + speedupText + '</td>';
                            html += '</tr>';
                        }
                    }

                    html += '</tbody></table>';
                }

                document.getElementById('results-tables').innerHTML = html;
            }

            // Create 2D plots and tables
            function create2DSection() {
                if (!ffta2dResults || !fftw2dResults) return;

                let html = '<h2 style="margin-top: 60px;">2D FFT Performance Visualizations</h2>';

                html += '<div class="plot-container">';
                html += '<h3>Runtime/Element vs Total Elements (All Categories)</h3>';
                html += '<div id="plot-2d-combined" class="plot"></div>';
                html += '</div>';

                html += '<div class="plot-container">';
                html += '<h3>Absolute Runtime (All Categories)</h3>';
                html += '<div id="plot-2d-absolute" class="plot"></div>';
                html += '</div>';

                html += '<h2>2D FFT Detailed Results</h2>';
                html += '<div id="results-tables-2d"></div>';

                document.getElementById('2d-section').innerHTML = html;

                // Create the plots
                create2DCombinedPlot();
                create2DAbsolutePlot();
                create2DResultsTables();
            }

            function create2DCombinedPlot() {
                if (!ffta2dResults || !fftw2dResults) return;

                const traces = [];

                for (const [catKey, catInfo] of Object.entries(categories2d)) {
                    const fftaCat = ffta2dResults.data.filter(d => d.category === catKey);
                    const fftwCat = fftw2dResults.data.filter(d => d.category === catKey);

                    if (fftaCat.length > 0) {
                        traces.push({
                            x: fftaCat.map(d => d.total_elements),
                            y: fftaCat.map(d => d.runtime_per_element * 1e9),
                            text: fftaCat.map(d => d.size_string),
                            name: 'FFTA: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, color: catInfo.color },
                            line: { width: 2, color: catInfo.color }
                        });
                    }

                    if (fftwCat.length > 0) {
                        traces.push({
                            x: fftwCat.map(d => d.total_elements),
                            y: fftwCat.map(d => d.runtime_per_element * 1e9),
                            text: fftwCat.map(d => d.size_string),
                            name: 'FFTW: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, symbol: 'square', color: catInfo.color },
                            line: { width: 2, dash: 'dash', color: catInfo.color }
                        });
                    }
                }

                const layout = {
                    xaxis: { title: 'Total Elements (rows × cols)', type: 'log' },
                    yaxis: { title: 'Runtime / Element (nanoseconds)', type: 'log' },
                    hovermode: 'closest',
                    showlegend: true,
                    legend: { x: 1.05, y: 1 }
                };

                Plotly.newPlot('plot-2d-combined', traces, layout, { responsive: true });
            }

            function create2DAbsolutePlot() {
                if (!ffta2dResults || !fftw2dResults) return;

                const traces = [];

                for (const [catKey, catInfo] of Object.entries(categories2d)) {
                    const fftaCat = ffta2dResults.data.filter(d => d.category === catKey);
                    const fftwCat = fftw2dResults.data.filter(d => d.category === catKey);

                    if (fftaCat.length > 0) {
                        traces.push({
                            x: fftaCat.map(d => d.total_elements),
                            y: fftaCat.map(d => d.median_time * 1e6),
                            text: fftaCat.map(d => d.size_string),
                            name: 'FFTA: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, color: catInfo.color },
                            line: { width: 2, color: catInfo.color }
                        });
                    }

                    if (fftwCat.length > 0) {
                        traces.push({
                            x: fftwCat.map(d => d.total_elements),
                            y: fftwCat.map(d => d.median_time * 1e6),
                            text: fftwCat.map(d => d.size_string),
                            name: 'FFTW: ' + catInfo.name,
                            type: 'scatter',
                            mode: 'lines+markers',
                            marker: { size: 8, symbol: 'square', color: catInfo.color },
                            line: { width: 2, dash: 'dash', color: catInfo.color }
                        });
                    }
                }

                const layout = {
                    xaxis: { title: 'Total Elements (rows × cols)', type: 'log' },
                    yaxis: { title: 'Runtime (microseconds)', type: 'log' },
                    hovermode: 'closest',
                    showlegend: true,
                    legend: { x: 1.05, y: 1 }
                };

                Plotly.newPlot('plot-2d-absolute', traces, layout, { responsive: true });
            }

            function create2DResultsTables() {
                if (!ffta2dResults || !fftw2dResults) return;

                let html = '';

                for (const [catKey, catInfo] of Object.entries(categories2d)) {
                    const fftaCat = ffta2dResults.data.filter(d => d.category === catKey);
                    const fftwCat = fftw2dResults.data.filter(d => d.category === catKey);

                    if (fftaCat.length === 0 && fftwCat.length === 0) continue;

                    html += '<h3>' + catInfo.name + '</h3>';
                    html += '<table>';
                    html += '<thead><tr>';
                    html += '<th>Array Size</th>';
                    html += '<th>Total Elements</th>';
                    html += '<th>FFTA Time (μs)</th>';
                    html += '<th>FFTW Time (μs)</th>';
                    html += '<th>FFTA Runtime/Element (ns)</th>';
                    html += '<th>FFTW Runtime/Element (ns)</th>';
                    html += '<th>Speedup</th>';
                    html += '</tr></thead>';
                    html += '<tbody>';

                    const allSizes = [...new Set([
                        ...fftaCat.map(d => d.size_string),
                        ...fftwCat.map(d => d.size_string)
                    ])].sort();

                    for (const sizeStr of allSizes) {
                        const fftaData = fftaCat.find(d => d.size_string === sizeStr);
                        const fftwData = fftwCat.find(d => d.size_string === sizeStr);

                        if (fftaData && fftwData) {
                            const fftaTime = (fftaData.median_time * 1e6).toFixed(3);
                            const fftwTime = (fftwData.median_time * 1e6).toFixed(3);
                            const fftaPerN = (fftaData.runtime_per_element * 1e9).toFixed(3);
                            const fftwPerN = (fftwData.runtime_per_element * 1e9).toFixed(3);
                            const speedup = fftwData.median_time / fftaData.median_time;
                            const speedupClass = speedup > 1 ? 'faster' : 'slower';
                            const speedupText = speedup > 1
                                ? speedup.toFixed(2) + 'x (FFTA faster)'
                                : (1/speedup).toFixed(2) + 'x (FFTW faster)';

                            html += '<tr>';
                            html += '<td>' + sizeStr + '</td>';
                            html += '<td>' + fftaData.total_elements + '</td>';
                            html += '<td>' + fftaTime + '</td>';
                            html += '<td>' + fftwTime + '</td>';
                            html += '<td>' + fftaPerN + '</td>';
                            html += '<td>' + fftwPerN + '</td>';
                            html += '<td class="' + speedupClass + '">' + speedupText + '</td>';
                            html += '</tr>';
                        }
                    }

                    html += '</tbody></table>';
                }

                document.getElementById('results-tables-2d').innerHTML = html;
            }

            // Initialize all plots on page load
            window.addEventListener('load', function() {
                createCombinedPlot();
                createAbsolutePlot();
                createCategoryPlot('odd_power_of_2', 'plot-odd');
                createCategoryPlot('even_power_of_2', 'plot-even');
                createCategoryPlot('power_of_3', 'plot-power3');
                createCategoryPlot('composite', 'plot-composite');
                createCategoryPlot('prime', 'plot-primes');
                createResultsTables();
                create2DSection();
            });
        </script>

        <footer style="margin-top: 40px; padding: 20px; text-align: center; color: #7f8c8d; border-top: 1px solid #ddd;">
            <p>FFTA.jl Performance Benchmark Suite</p>
            <p>Interactive plots powered by Plotly.js</p>
        </footer>
    </body>
    </html>
    """

    # Write HTML file
    output_file = joinpath(@__DIR__, "benchmark_report.html")
    write(output_file, html)
    println("HTML report saved to: ", output_file)
    println("Open in a web browser to view interactive Plotly.js charts")

    return output_file
end

# Generate report
generate_html_report()
