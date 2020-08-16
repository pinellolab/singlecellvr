AFRAME.registerComponent('loading', {
    schema: {
        time: {type: 'number', default: 10},
        show: {type: 'boolean', default: true},
        zDepth: {type: 'number', default: -1}
    },
    init: function() {
        this.initialized = performance.now();
        const width = this.visibleWidthAtZDepth(this.data.zDepth);
        const height = this.visibleHeightAtZDepth(this.data.zDepth);
        this.el.setAttribute('width', width);
        this.el.setAttribute('height', height);

    },
    update: function(oldData) {
        if (this.data.show !== oldData.show && !this.data.show) {
            const dismissTime = performance.now();
            const elapsedSeconds = (dismissTime - this.initialized) / 1000;
            if (elapsedSeconds > this.data.time) {
                this.el.object3D.visible = false;
                document.getElementById("hud").setAttribute('visible', true);
                document.getElementById("player-camera").setAttribute("look-controls", "");
            } else {
                setTimeout(() => {
                    this.el.object3D.visible = false;
                    document.getElementById("hud").setAttribute('visible', true);
                    document.getElementById("player-camera").setAttribute("look-controls", "");
                }, ((this.data.time - elapsedSeconds) * 1000))
            }
        }
    },
    tick: function() {
        const cameraPosition = document.getElementById('rig').object3D.position;
        this.el.object3D.position.set(cameraPosition.x, cameraPosition.y, cameraPosition.z - 1);
    },
    // These are duplicated, fix later
    visibleHeightAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        // compensate for cameras not positioned at z=0
        const cameraOffset = camera.position.z;
        if ( depth < cameraOffset ) depth -= cameraOffset;
        else depth += cameraOffset;
      
        // vertical fov in radians
        const vFOV = camera.fov * Math.PI / 180;
      
        // Math.abs to ensure the result is always positive
        return 2 * Math.tan( vFOV / 2 ) * Math.abs( depth );
    },
    visibleWidthAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        const height = this.visibleHeightAtZDepth( depth, camera );
        let width = height * camera.aspect;
        return width;
    }

});