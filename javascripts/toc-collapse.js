/**
 * Collapsible Table of Contents (TOC) functionality
 * Makes top-level TOC sections collapsible on click
 */

document.addEventListener("DOMContentLoaded", function() {
    // Wait for the page to fully load
    setTimeout(function() {
        initCollapsibleTOC();
    }, 100);
});

function initCollapsibleTOC() {
    // Find the TOC navigation in the right sidebar
    const tocNav = document.querySelector('.md-sidebar--secondary .md-nav--secondary');
    
    if (!tocNav) return;
    
    // Find all nested items (items with sub-items)
    const nestedItems = tocNav.querySelectorAll('.md-nav__item--nested');
    
    nestedItems.forEach(function(item) {
        const link = item.querySelector(':scope > .md-nav__link');
        
        if (!link) return;
        
        // Add click handler to toggle collapse
        link.addEventListener('click', function(e) {
            // Don't prevent default if it's an actual link with href
            if (this.getAttribute('href') && this.getAttribute('href') !== '#') {
                return;
            }
            
            e.preventDefault();
            e.stopPropagation();
            
            // Toggle collapsed state
            item.classList.toggle('collapsed');
        });
        
        // Optional: Start with some sections collapsed
        // Uncomment the next line to start all sections collapsed
        // item.classList.add('collapsed');
    });
}

// Re-initialize on page navigation (for instant loading)
if (typeof document$ !== 'undefined') {
    document$.subscribe(function() {
        initCollapsibleTOC();
    });
}
