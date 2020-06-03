/**
   * Remove stubborn lighting from environment
   */
  AFRAME.registerComponent('remove-lighting', {
    
    dependencies: ['environment'],
    
    update: function() {
      var lights = this.el.querySelectorAll('[light]');
      
      for (var i = 0; i < lights.length; i++) {
        lights[i].removeAttribute('light');
      }
    }
    
});